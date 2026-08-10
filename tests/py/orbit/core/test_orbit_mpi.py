import pytest
from orbit.core.orbit_mpi import (
    mpi_comm,
    mpi_datatype,
    mpi_op,
    MPI_Allreduce,
    MPI_Comm_size,
)

_MPI_SIZE = MPI_Comm_size(mpi_comm.MPI_COMM_WORLD)

_OPS = [
    pytest.param(mpi_op.MPI_SUM, "sum", id="MPI_SUM"),
    pytest.param(mpi_op.MPI_MAX, "max", id="MPI_MAX"),
    pytest.param(mpi_op.MPI_MIN, "min", id="MPI_MIN"),
    pytest.param(mpi_op.MPI_PROD, "prod", id="MPI_PROD"),
]


def _expected(val, op_name, size):
    if op_name == "sum":
        return val * size
    elif op_name == "max":
        return val
    elif op_name == "min":
        return val
    elif op_name == "prod":
        return val**size
    raise ValueError(f"unknown op: {op_name}")


def _expected_seq(vals, op_name, size):
    return tuple(_expected(v, op_name, size) for v in vals)


class TestMPI_Allreduce:
    @pytest.mark.parametrize("op,op_name", _OPS)
    def test_double_scalar(self, op, op_name):
        val = 3.14
        expected = _expected(val, op_name, _MPI_SIZE)
        res = MPI_Allreduce(val, mpi_datatype.MPI_DOUBLE, op, mpi_comm.MPI_COMM_WORLD)
        assert isinstance(res, float)
        assert res == pytest.approx(expected)

    @pytest.mark.parametrize("op,op_name", _OPS)
    def test_int_scalar(self, op, op_name):
        val = 42
        expected = _expected(val, op_name, _MPI_SIZE)
        res = MPI_Allreduce(val, mpi_datatype.MPI_INT, op, mpi_comm.MPI_COMM_WORLD)
        assert isinstance(res, int)
        assert res == expected

    @pytest.mark.parametrize("op,op_name", _OPS)
    def test_double_tuple(self, op, op_name):
        vals = (1.5, 2.5, 3.5)
        expected = _expected_seq(vals, op_name, _MPI_SIZE)
        res = MPI_Allreduce(vals, mpi_datatype.MPI_DOUBLE, op, mpi_comm.MPI_COMM_WORLD)
        assert isinstance(res, tuple)
        assert len(res) == len(expected)
        for a, b in zip(res, expected):
            assert a == pytest.approx(b)

    @pytest.mark.parametrize("op,op_name", _OPS)
    def test_int_tuple(self, op, op_name):
        vals = (10, 20, 30)
        expected = _expected_seq(vals, op_name, _MPI_SIZE)
        res = MPI_Allreduce(vals, mpi_datatype.MPI_INT, op, mpi_comm.MPI_COMM_WORLD)
        assert isinstance(res, tuple)
        assert list(res) == list(expected)

    @pytest.mark.parametrize("op,op_name", _OPS)
    def test_double_list(self, op, op_name):
        vals = [1.5, 2.5, 3.5]
        expected = _expected_seq(vals, op_name, _MPI_SIZE)
        res = MPI_Allreduce(vals, mpi_datatype.MPI_DOUBLE, op, mpi_comm.MPI_COMM_WORLD)
        assert isinstance(res, tuple)
        assert len(res) == len(expected)
        for a, b in zip(res, expected):
            assert a == pytest.approx(b)

    @pytest.mark.parametrize("op,op_name", _OPS)
    def test_int_list(self, op, op_name):
        vals = [10, 20, 30]
        expected = _expected_seq(vals, op_name, _MPI_SIZE)
        res = MPI_Allreduce(vals, mpi_datatype.MPI_INT, op, mpi_comm.MPI_COMM_WORLD)
        assert isinstance(res, tuple)
        assert list(res) == list(expected)

    @pytest.mark.skip(
        reason="numpy arrays not yet supported by MPI_Allreduce C wrapper – segfaults"
    )
    @pytest.mark.parametrize("op,op_name", _OPS)
    def test_double_ndarray(self, op, op_name):
        np = pytest.importorskip("numpy")
        vals = np.array([1.5, 2.5, 3.5])
        expected = _expected_seq(vals, op_name, _MPI_SIZE)
        res = MPI_Allreduce(vals, mpi_datatype.MPI_DOUBLE, op, mpi_comm.MPI_COMM_WORLD)
        assert isinstance(res, tuple)
        assert len(res) == len(expected)
        for a, b in zip(res, expected):
            assert a == pytest.approx(float(b))

    @pytest.mark.skip(
        reason="numpy arrays not yet supported by MPI_Allreduce C wrapper – segfaults"
    )
    @pytest.mark.parametrize("op,op_name", _OPS)
    def test_int_ndarray(self, op, op_name):
        np = pytest.importorskip("numpy")
        vals = np.array([10, 20, 30], dtype=np.int64)
        expected = _expected_seq(vals, op_name, _MPI_SIZE)
        res = MPI_Allreduce(vals, mpi_datatype.MPI_INT, op, mpi_comm.MPI_COMM_WORLD)
        assert isinstance(res, tuple)
        assert list(res) == list(map(int, expected))
