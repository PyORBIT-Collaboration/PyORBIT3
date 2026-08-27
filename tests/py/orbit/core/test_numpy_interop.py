import os
import shutil
import subprocess
import sys
import textwrap

import numpy as np
import pytest

from orbit.core.bunch import Bunch
from orbit.core.spacecharge import Grid3D


NUMPY_INTEROP_ENABLED = hasattr(Bunch, "to_numpy") and hasattr(
    Grid3D, "to_numpy"
)

pytestmark = pytest.mark.skipif(
    not NUMPY_INTEROP_ENABLED,
    reason="PyORBIT was built without experimental NumPy interoperability",
)


def _process_failure(result):
    return textwrap.dedent(
        f"""
        subprocess exited with status {result.returncode}
        stdout:
        {result.stdout}
        stderr:
        {result.stderr}
        """
    )


def _assert_missing_numpy_raises_import_error(extension_name):
    """NumPy C-API initialization failure must not terminate Python."""
    script = textwrap.dedent(
        f"""
        import sys

        class BlockNumpyImports:
            @staticmethod
            def find_spec(fullname, path=None, target=None):
                if fullname == "numpy" or fullname.startswith("numpy."):
                    raise ModuleNotFoundError("NumPy is unavailable for this test")
                return None

        sys.meta_path.insert(0, BlockNumpyImports())
        try:
            from orbit.core import {extension_name}
        except ImportError:
            raise SystemExit(0)
        raise SystemExit("expected importing without NumPy to fail")
        """
    )
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)

    result = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
        env=env,
    )

    assert result.returncode == 0, _process_failure(result)


class TestBunchNumpyInterop:
    def test_update_from_numpy_is_atomic_on_invalid_input(self):
        """A rejected update must leave the original particles untouched."""
        original = np.arange(18, dtype=np.float64).reshape(3, 6)
        bunch = Bunch.from_numpy(original)

        with pytest.raises(ValueError, match=r"shape \(nparts, 6\)"):
            bunch.update_from_numpy(np.zeros((4, 5), dtype=np.float64))

        np.testing.assert_array_equal(bunch.to_numpy(), original)

    def test_to_numpy_gather_on_single_rank(self):
        """Gathering on a serial communicator returns the local array."""
        expected = np.arange(18, dtype=np.float64).reshape(3, 6)
        bunch = Bunch.from_numpy(expected)

        np.testing.assert_array_equal(bunch.to_numpy(gather=True), expected)

    def test_missing_numpy_raises_import_error_instead_of_aborting(self):
        _assert_missing_numpy_raises_import_error("bunch")

    def test_update_from_numpy_uses_only_local_particle_count(self, tmp_path):
        """Each MPI rank must clear only its locally allocated particles."""
        mpirun = shutil.which("mpirun") or shutil.which("mpiexec")
        if mpirun is None:
            pytest.skip("MPI launcher is not available")

        script_path = tmp_path / "numpy_bunch_mpi_update.py"
        script_path.write_text(
            textwrap.dedent(
                """
                import numpy as np

                from orbit.core.bunch import Bunch
                from orbit.core.orbit_mpi import (
                    MPI_Comm_rank,
                    MPI_Comm_size,
                    mpi_comm,
                )

                comm = mpi_comm.MPI_COMM_WORLD
                rank = MPI_Comm_rank(comm)
                size = MPI_Comm_size(comm)
                if size != 2:
                    print("PYORBIT_MPI_DISABLED", flush=True)
                    raise SystemExit(0)

                # Rank zero deliberately has only the minimum allocation while
                # rank one makes the global count much larger than that allocation.
                local_size = 0 if rank == 0 else 100_000
                bunch = Bunch.from_numpy(
                    np.ones((local_size, 6), dtype=np.float64)
                )
                bunch.update_from_numpy(
                    np.full((1, 6), rank, dtype=np.float64)
                )

                assert bunch.getSize() == 1
                np.testing.assert_array_equal(
                    bunch.to_numpy(), np.full((1, 6), rank, dtype=np.float64)
                )
                """
            )
        )

        env = os.environ.copy()
        env.setdefault("OMPI_ALLOW_RUN_AS_ROOT", "1")
        env.setdefault("OMPI_ALLOW_RUN_AS_ROOT_CONFIRM", "1")
        result = subprocess.run(
            [mpirun, "-np", "2", sys.executable, str(script_path)],
            capture_output=True,
            text=True,
            timeout=60,
            check=False,
            env=env,
        )

        if "PYORBIT_MPI_DISABLED" in result.stdout:
            pytest.skip("PyORBIT was built without MPI support")

        assert result.returncode == 0, _process_failure(result)

    def test_to_numpy_gathers_particles_across_mpi_ranks(self, tmp_path):
        """Gathered rows are rank ordered and returned only on the chosen root."""
        mpirun = shutil.which("mpirun") or shutil.which("mpiexec")
        if mpirun is None:
            pytest.skip("MPI launcher is not available")

        script_path = tmp_path / "numpy_bunch_mpi_gather.py"
        script_path.write_text(
            textwrap.dedent(
                """
                import numpy as np
                import pytest

                from orbit.core.bunch import Bunch
                from orbit.core.orbit_mpi import (
                    MPI_Comm_rank,
                    MPI_Comm_size,
                    mpi_comm,
                )

                comm = mpi_comm.MPI_COMM_WORLD
                rank = MPI_Comm_rank(comm)
                size = MPI_Comm_size(comm)
                if size != 3:
                    print("PYORBIT_MPI_DISABLED", flush=True)
                    raise SystemExit(0)

                local_sizes = (2, 0, 3)

                def particles_for_rank(mpi_rank):
                    count = local_sizes[mpi_rank]
                    return (
                        np.arange(count * 6, dtype=np.float64).reshape(count, 6)
                        + 100.0 * mpi_rank
                    )

                local_particles = particles_for_rank(rank)
                expected_global = np.concatenate(
                    [particles_for_rank(mpi_rank) for mpi_rank in range(size)]
                )
                bunch = Bunch.from_numpy(local_particles)

                # Complete every collective before making assertions so an
                # assertion failure on one rank cannot strand another rank.
                default_root_result = bunch.to_numpy(gather=True)
                empty_root_result = bunch.to_numpy(gather=True, root=1)

                empty_bunch = Bunch.from_numpy(
                    np.empty((0, 6), dtype=np.float64)
                )
                empty_global_result = empty_bunch.to_numpy(gather=True, root=2)

                # Gathering must not change the ordinary local conversion.
                np.testing.assert_array_equal(bunch.to_numpy(), local_particles)

                if rank == 0:
                    np.testing.assert_array_equal(
                        default_root_result, expected_global
                    )
                else:
                    assert default_root_result is None

                # Rank one owns no local particles but can still be the root.
                if rank == 1:
                    np.testing.assert_array_equal(
                        empty_root_result, expected_global
                    )
                else:
                    assert empty_root_result is None

                if rank == 2:
                    assert empty_global_result.shape == (0, 6)
                    assert empty_global_result.dtype == np.float64
                else:
                    assert empty_global_result is None

                with pytest.raises(
                    ValueError, match="root must be between 0 and 2"
                ):
                    bunch.to_numpy(gather=True, root=3)
                """
            )
        )

        env = os.environ.copy()
        env.setdefault("OMPI_ALLOW_RUN_AS_ROOT", "1")
        env.setdefault("OMPI_ALLOW_RUN_AS_ROOT_CONFIRM", "1")
        result = subprocess.run(
            [mpirun, "-np", "3", sys.executable, str(script_path)],
            capture_output=True,
            text=True,
            timeout=60,
            check=False,
            env=env,
        )

        if "PYORBIT_MPI_DISABLED" in result.stdout:
            pytest.skip("PyORBIT was built without MPI support")

        assert result.returncode == 0, _process_failure(result)


class TestGrid3DNumpyInterop:
    def test_to_numpy_round_trips_rectangular_grid(self):
        """Exporting a non-square grid must not access memory out of bounds."""
        script = textwrap.dedent(
            """
            import numpy as np

            from orbit.core.spacecharge import Grid3D

            expected = np.arange(
                4 * 5 * 6, dtype=np.float64
            ).reshape(4, 5, 6)
            grid = Grid3D(4, 5, 6)
            grid.from_numpy(expected)
            np.testing.assert_array_equal(grid.to_numpy(), expected)
            """
        )

        result = subprocess.run(
            [sys.executable, "-c", script],
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )

        assert result.returncode == 0, _process_failure(result)

    def test_from_numpy_uses_xyz_order(self):
        values_xyz = np.arange(4 * 5 * 6, dtype=np.float64).reshape(4, 5, 6)
        grid = Grid3D(4, 5, 6)

        grid.from_numpy(values_xyz)

        for ix, iy, iz in ((0, 0, 0), (1, 2, 3), (3, 4, 5)):
            assert grid.getValueOnGrid(ix, iy, iz) == values_xyz[ix, iy, iz]

    def test_missing_numpy_raises_import_error_instead_of_aborting(self):
        _assert_missing_numpy_raises_import_error("spacecharge")
