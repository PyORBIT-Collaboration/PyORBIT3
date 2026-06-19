python = import('python').find_installation('python3', pure: false)

cpp = meson.get_compiler('cpp')

if cpp.has_argument('-march=native')
  add_project_arguments('-march=native', language: 'cpp')
endif

dependencies = []

dependencies += python.dependency()
dependencies += dependency('fftw3', version: '>= 3.0.0', required: true)

py_link_args = []

if host_machine.system() == 'darwin'
  # The unresolved symbols will be resolved once the Python interpreter loads libpyorbit.
  # The linker on MacOS is a bit more strict, so we need to explicitly tell it that
  # symbols from libpython will be nont resolvable for now.
  py_link_args += ['-Wl,-undefined,dynamic_lookup']
endif

# Detecting if MPICH or OPENMPI are installed and enabling support if present
mpi_use = get_option('USE_MPI')
use_mpi_flag = '0'

if mpi_use == 'mpich'
  message('Requested to use MPICH as the MPI implementation.')
  dependencies += dependency('mpich', version: '>= 4.0.0', required: true)
  add_project_arguments('-DUSE_MPI=1', language: 'cpp')
  use_mpi_flag = '1'
elif mpi_use == 'ompi'
  message('Requested to use OpenMPI as the MPI implementation.')
  dependencies += dependency('ompi', version: '>= 4.0.0', required: true)
  add_project_arguments('-DUSE_MPI=1', language: 'cpp')
  use_mpi_flag = '1'
elif mpi_use == 'none'
  message('Requested to not use MPI.')
else
   mpich_dependency = dependency('mpich', version: '>= 4.0.0', required: false)
   openmpi_dependency = dependency('ompi', version: '>= 4.0.0', required: false)

   if mpich_dependency.found()
     add_project_arguments('-DUSE_MPI=1', language: 'cpp')
     dependencies += mpich_dependency
     message('Using MPICH as the MPI implementation.')
     use_mpi_flag = '1'
   elif openmpi_dependency.found()
     add_project_arguments('-DUSE_MPI=1', language: 'cpp')
     dependencies += openmpi_dependency
     message('Using OpenMPI as the MPI implementation.')
     use_mpi_flag = '1'
   else
     message('MPI will not be used.')
   endif

endif

sources = files([
  'linac/tracking/linac_tracking.cc',
  'linac/rfgap/BaseRfGap_slow.cc',
  'linac/rfgap/BaseRfGap.cc',
  'linac/rfgap/RfGapThreePointTTF.cc',
  'linac/rfgap/MatrixRfGap.cc',
  'linac/rfgap/RfGapTTF.cc',
  'linac/rfgap/SuperFishFieldSource.cc',
  'linac/rfgap/RfGapTTF_slow.cc',
  'linac/rfgap/RfGapThreePointTTF_slow.cc',
  'orbit/OrbitConst.cc',
  'orbit/Bunch.cc',
  'orbit/SyncPart.cc',
  'orbit/BunchDiagnostics/BunchTwissAnalysis.cc',
  'orbit/BunchDiagnostics/BunchTuneAnalysis.cc',
  'orbit/FieldTracker/FieldTracker.cc',
  'orbit/MaterialInteractions/Foil.cc',
  'orbit/MaterialInteractions/Collimator.cc',
  'orbit/MaterialInteractions/numrecipes.cc',
  'orbit/MaterialInteractions/cross_sections.cc',
  'orbit/MaterialInteractions/MaterialInteractions.cc',
  'orbit/RFCavities/Dual_Harmonic_Cav.cc',
  'orbit/RFCavities/Frequency_Cav.cc',
  'orbit/RFCavities/Barrier_Cav.cc',
  'orbit/RFCavities/Harmonic_Cav.cc',
  'orbit/Impedances/TImpedance.cc',
  'orbit/Impedances/LImpedance.cc',
  'orbit/SynchPartRedefinition/SynchPartRedefinitionZdE.cc',
  'orbit/Errors/errorbase.cc',
  'orbit/ParticlesAttributes/ParticleInitialCoordinates.cc',
  'orbit/ParticlesAttributes/ParticleIdNumber.cc',
  'orbit/ParticlesAttributes/ParticleMacroSize.cc',
  'orbit/ParticlesAttributes/ParticleAttributesFactory.cc',
  'orbit/ParticlesAttributes/WaveFunctionAmplitudes.cc',
  'orbit/ParticlesAttributes/pq_coordinates.cc',
  'orbit/ParticlesAttributes/Evolution.cc',
  'orbit/ParticlesAttributes/ParticleAttributes.cc',
  'orbit/ParticlesAttributes/TurnNumberAttributes.cc',
  'orbit/ParticlesAttributes/LostParticleAttributes.cc',
  'orbit/ParticlesAttributes/ParticlePhaseAttributes.cc',
  'orbit/ParticlesAttributes/part_time.cc',
  'orbit/ParticlesAttributes/AtomPopulations.cc',
  'orbit/Apertures/CircleApertureShape.cc',
  'orbit/Apertures/PhaseAperture.cc',
  'orbit/Apertures/Aperture.cc',
  'orbit/Apertures/RectangularApertureShape.cc',
  'orbit/Apertures/PyBaseApertureShape.cc',
  'orbit/Apertures/ConvexApertureShape.cc',
  'orbit/Apertures/EnergyAperture.cc',
  'orbit/Apertures/BaseAperture.cc',
  'orbit/Apertures/EllipseApertureShape.cc',
  'orbit/Apertures/CompositeApertureShape.cc',
  'orbit/Apertures/BaseApertureShape.cc',
  'trackerrk4/ExternalEffects.cc',
  'trackerrk4/ExtEffectsContainer.cc',
  'trackerrk4/PyExternalEffects.cc',
  'trackerrk4/RungeKuttaTracker.cc',
  'mpi/orbit_mpi.cc',
  'spacecharge/PoissonSolverFFT3D.cc',
  'spacecharge/SpaceChargeCalcUnifEllipse.cc',
  'spacecharge/ForceSolver2D.cc',
  'spacecharge/PoissonSolverFFT2D.cc',
  'spacecharge/PoissonSolver2D.cc',
  'spacecharge/SpaceChargeForceCalc2p5D.cc',
  'spacecharge/PoissonSolver3D.cc',
  'spacecharge/SpaceChargeCalc3D.cc',
  'spacecharge/ForceSolverFFT2D.cc',
  'spacecharge/SpaceChargeCalcSliceBySlice2D.cc',
  'spacecharge/Grid3D.cc',
  'spacecharge/SpaceChargeCalc2p5Drb.cc',
  'spacecharge/Grid2D.cc',
  'spacecharge/SpaceChargeCalc2p5D.cc',
  'spacecharge/Grid1D.cc',
  'spacecharge/LSpaceChargeCalc.cc',
  'spacecharge/ShapedBoundary2D.cc',
  'spacecharge/BaseBoundary2D.cc',
  'spacecharge/UniformEllipsoidFieldCalculator.cc',
  'utils/FieldSourceContainer.cc',
  'utils/AttributesBucket.cc',
  'utils/bessel.cc',
  'utils/PyBaseFieldSource.cc',
  'utils/OU_SplineCH.cc',
  'utils/CppPyWrapper.cc',
  'utils/StringUtils.cc',
  'utils/Random.cc',
  'utils/BaseFieldSource.cc',
  'utils/OU_Function.cc',
  'utils/MathPolynomial.cc',
  'utils/BufferStore.cc',
  'utils/statistics/StatMoments2D.cc',
  'utils/field_sources/LoopFieldSource.cc',
  'utils/field_sources/QuadFieldSource.cc',
  'utils/field_sources/ShiftedFieldSource.cc',
  'utils/field_sources/DipoleFieldSource.cc',
  'utils/field_sources/MagnetFieldSourceGrid3D.cc',
  'utils/integration/gauss_legendre_points.cc',
  'utils/integration/GaussLegendreIntegrator.cc',
  'utils/bunch/ParticlesWithIdFunctions.cc',
  'utils/bunch/TwissFilteringFunctions.cc',
  'utils/bunch/InitialCoordsAttrFunctions.cc',
  'utils/bunch/BunchExtremaCalculator.cc',
  'utils/harmonic_analysis/HarmonicData.cc',
  'utils/polynomial/OU_Polynomial.cc',
  'utils/matrix/PhaseVector.cc',
  'utils/matrix/MatrixOperations.cc',
  'utils/matrix/Matrix.cc',
  'teapot/teapotbase.cc',
  'teapot/MatrixGenerator.cc'
])

inc = include_directories([
  'main',
  'mpi',
  'orbit',
  'orbit/Apertures',
  'orbit/BunchDiagnostics',
  'orbit/Errors',
  'orbit/FieldTracker',
  'orbit/Impedances',
  'orbit/MaterialInteractions',
  'orbit/ParticlesAttributes',
  'orbit/RFCavities',
  'orbit/SynchPartRedefinition',
  'spacecharge',
  'teapot',
  'trackerrk4',
  'utils',
  'utils/bunch',
  'utils/ellipticalint',
  'utils/field_sources',
  'utils/harmonic_analysis',
  'utils/integration',
  'utils/matrix',
  'utils/polynomial',
  'utils/statistics',
  'linac',
  'linac/rfgap',
  'linac/tracking',
])

core_lib = library(
  'pyorbit3',
  sources: sources,
  include_directories: inc,
  dependencies: dependencies,
  link_args: py_link_args,
  install: true,
)

core_dep = declare_dependency(link_with : core_lib, dependencies: dependencies)

api_headers = [
  'linac/rfgap/BaseRfGap.hh',
  'linac/rfgap/BaseRfGap_slow.hh',
  'linac/rfgap/MatrixRfGap.hh',
  'linac/rfgap/RfGapThreePointTTF.hh',
  'linac/rfgap/RfGapThreePointTTF_slow.hh',
  'linac/rfgap/RfGapTTF.hh',
  'linac/rfgap/RfGapTTF_slow.hh',
  'linac/rfgap/SuperFishFieldSource.hh',
  'linac/tracking/linac_tracking.hh',
  'main/pyORBIT_Object.hh',
  'mpi/orbit_mpi.hh',
  'orbit/Apertures/Aperture.hh',
  'orbit/Apertures/BaseAperture.hh',
  'orbit/Apertures/BaseApertureShape.hh',
  'orbit/Apertures/CircleApertureShape.hh',
  'orbit/Apertures/CompositeApertureShape.hh',
  'orbit/Apertures/ConvexApertureShape.hh',
  'orbit/Apertures/EllipseApertureShape.hh',
  'orbit/Apertures/EnergyAperture.hh',
  'orbit/Apertures/PhaseAperture.hh',
  'orbit/Apertures/PyBaseApertureShape.hh',
  'orbit/Apertures/RectangularApertureShape.hh',
  'orbit/Bunch.hh',
  'orbit/BunchDiagnostics/BunchTuneAnalysis.hh',
  'orbit/BunchDiagnostics/BunchTwissAnalysis.hh',
  'orbit/Errors/errorbase.hh',
  'orbit/FieldTracker/FieldTracker.hh',
  'orbit/Impedances/LImpedance.hh',
  'orbit/Impedances/TImpedance.hh',
  'orbit/MaterialInteractions/Collimator.hh',
  'orbit/MaterialInteractions/Foil.hh',
  'orbit/MaterialInteractions/MaterialInteractions.hh',
  'orbit/MaterialInteractions/cross_sections.hh',
  'orbit/MaterialInteractions/numrecipes.hh',
  'orbit/OrbitConst.hh',
  'orbit/ParticlesAttributes/AtomPopulations.hh',
  'orbit/ParticlesAttributes/Evolution.hh',
  'orbit/ParticlesAttributes/LostParticleAttributes.hh',
  'orbit/ParticlesAttributes/ParticleAttributes.hh',
  'orbit/ParticlesAttributes/ParticleAttributesFactory.hh',
  'orbit/ParticlesAttributes/ParticleIdNumber.hh',
  'orbit/ParticlesAttributes/ParticleInitialCoordinates.hh',
  'orbit/ParticlesAttributes/ParticleMacroSize.hh',
  'orbit/ParticlesAttributes/ParticlePhaseAttributes.hh',
  'orbit/ParticlesAttributes/TurnNumberAttributes.hh',
  'orbit/ParticlesAttributes/WaveFunctionAmplitudes.hh',
  'orbit/ParticlesAttributes/part_time.hh',
  'orbit/ParticlesAttributes/pq_coordinates.hh',
  'orbit/RFCavities/Barrier_Cav.hh',
  'orbit/RFCavities/Dual_Harmonic_Cav.hh',
  'orbit/RFCavities/Frequency_Cav.hh',
  'orbit/RFCavities/Harmonic_Cav.hh',
  'orbit/SynchPartRedefinition/SynchPartRedefinitionZdE.hh',
  'orbit/SyncPart.hh',
  'spacecharge/BaseBoundary2D.hh',
  'spacecharge/ForceSolver2D.hh',
  'spacecharge/ForceSolverFFT2D.hh',
  'spacecharge/Grid1D.hh',
  'spacecharge/Grid2D.hh',
  'spacecharge/Grid3D.hh',
  'spacecharge/LSpaceChargeCalc.hh',
  'spacecharge/PoissonSolver2D.hh',
  'spacecharge/PoissonSolver3D.hh',
  'spacecharge/PoissonSolverFFT2D.hh',
  'spacecharge/PoissonSolverFFT3D.hh',
  'spacecharge/ShapedBoundary2D.hh',
  'spacecharge/SpaceChargeCalc2p5D.hh',
  'spacecharge/SpaceChargeCalc2p5Drb.hh',
  'spacecharge/SpaceChargeCalc3D.hh',
  'spacecharge/SpaceChargeCalcSliceBySlice2D.hh',
  'spacecharge/SpaceChargeCalcUnifEllipse.hh',
  'spacecharge/SpaceChargeForceCalc2p5D.hh',
  'spacecharge/UniformEllipsoidFieldCalculator.hh',
  'teapot/MatrixGenerator.hh',
  'teapot/teapotbase.hh',
  'trackerrk4/ExtEffectsContainer.hh',
  'trackerrk4/ExternalEffects.hh',
  'trackerrk4/PyExternalEffects.hh',
  'trackerrk4/RungeKuttaTracker.hh',
  'utils/AttributesBucket.hh',
  'utils/BaseFieldSource.hh',
  'utils/BufferStore.hh',
  'utils/CppPyWrapper.hh',
  'utils/FieldSourceContainer.hh',
  'utils/MathPolynomial.hh',
  'utils/OU_Function.hh',
  'utils/OU_SplineCH.hh',
  'utils/PyBaseFieldSource.hh',
  'utils/Random.hh',
  'utils/StringUtils.hh',
  'utils/bessel.hh',
  'utils/bunch/BunchExtremaCalculator.hh',
  'utils/bunch/InitialCoordsAttrFunctions.hh',
  'utils/bunch/ParticlesWithIdFunctions.hh',
  'utils/bunch/TwissFilteringFunctions.hh',
  'utils/ellipticalint/elliptint.hh',
  'utils/field_sources/DipoleFieldSource.hh',
  'utils/field_sources/LoopFieldSource.hh',
  'utils/field_sources/MagnetFieldSourceGrid3D.hh',
  'utils/field_sources/QuadFieldSource.hh',
  'utils/field_sources/ShiftedFieldSource.hh',
  'utils/harmonic_analysis/HarmonicData.hh',
  'utils/integration/GaussLegendreIntegrator.hh',
  'utils/integration/gauss_legendre_points.hh',
  'utils/matrix/Matrix.hh',
  'utils/matrix/MatrixOperations.hh',
  'utils/matrix/PhaseVector.hh',
  'utils/polynomial/OU_Polynomial.hh',
  'utils/statistics/StatMoments2D.hh',
  'utils/tcomplex.hh',
]

install_headers(api_headers, subdir : 'pyorbit3', preserve_path : true)

pc_requires = ['fftw3 >= 3.0.0']
if get_option('USE_MPI') == 'mpich'
  pc_requires += ['mpich >= 4.0.0']
elif get_option('USE_MPI') == 'ompi'
  pc_requires += ['ompi >= 4.0.0']
endif

pkgconfig = import('pkgconfig')

pkgconfig.generate(
  name: 'pyorbit3',
  description: 'PyORBIT3 core library',
  filebase: 'pyorbit3',
  version: meson.project_version(),
  libraries: core_lib,
  requires: pc_requires,
  subdirs: 'pyorbit3',
  extra_cflags: '-DUSE_MPI=' + use_mpi_flag,
  install_dir: get_option('libdir') / 'pkgconfig',
)

ext_modules = {
  'orbit_mpi': files(
    'core/mpi_init.cc',
    'mpi/wrap_orbit_mpi.cc',
    'mpi/wrap_mpi_comm.cc',
    'mpi/wrap_mpi_datatype.cc',
    'mpi/wrap_mpi_request.cc',
    'mpi/wrap_mpi_op.cc',
    'mpi/wrap_mpi_status.cc',
    'mpi/wrap_mpi_group.cc',
  ),
  'bunch': files(
    'core/bunch_init.cc',
    'core/type_accessors.cc',
    'orbit/wrap_bunch.cc',
    'orbit/wrap_syncpart.cc',
    'orbit/BunchDiagnostics/wrap_bunch_tune_analysis.cc',
    'orbit/BunchDiagnostics/wrap_bunch_twiss_analysis.cc',
    'orbit/SynchPartRedefinition/wrap_synch_part_redefinition_z_de.cc',
  ),
  'spacecharge': files(
    'core/spacecharge_init.cc',
    'core/type_accessors.cc',
    'spacecharge/wrap_spacecharge.cc',
    'spacecharge/wrap_grid1D.cc',
    'spacecharge/wrap_grid2D.cc',
    'spacecharge/wrap_grid3D.cc',
    'spacecharge/wrap_poissonsolverfft2d.cc',
    'spacecharge/wrap_poissonsolverfft3d.cc',
    'spacecharge/wrap_forcesolverfft2d.cc',
    'spacecharge/wrap_boundary2d.cc',
    'spacecharge/wrap_spacechargecalc2p5d.cc',
    'spacecharge/wrap_spacechargeforcecalc2p5d.cc',
    'spacecharge/wrap_spacechargecalc2p5d_rb.cc',
    'spacecharge/wrap_spacechargecalc_slicebyslice_2D.cc',
    'spacecharge/wrap_lspacechargecalc.cc',
    'spacecharge/wrap_spacechargecalc3d.cc',
    'spacecharge/wrap_uniform_ellipsoid_field_calculator.cc',
    'spacecharge/wrap_spacechargecalc_uniform_ellipse.cc',
  ),
  'trackerrk4': files(
    'core/trackerrk4_init.cc',
    'trackerrk4/wrap_trackerrk4.cc',
    'trackerrk4/wrap_runge_kutta_tracker.cc',
    'trackerrk4/wrap_py_external_effects.cc',
    'trackerrk4/wrap_ext_effects_container.cc',
  ),
  'teapot_base': files(
    'core/teapot_base_init.cc',
    'core/type_accessors.cc',
    'teapot/wrap_teapotbase.cc',
    'teapot/wrap_matrix_generator.cc',
  ),
  'linac': files(
    'core/linac_init.cc',
    'core/type_accessors.cc',
    'linac/wrap_linacmodule.cc',
    'linac/tracking/wrap_linac_tracking.cc',
    'linac/rfgap/wrap_BaseRfGap.cc',
    'linac/rfgap/wrap_BaseRfGap_slow.cc',
    'linac/rfgap/wrap_MatrixRfGap.cc',
    'linac/rfgap/wrap_RfGapTTF.cc',
    'linac/rfgap/wrap_RfGapTTF_slow.cc',
    'linac/rfgap/wrap_SuperFishFieldSource.cc',
    'linac/rfgap/wrap_RfGapThreePointTTF.cc',
    'linac/rfgap/wrap_RfGapThreePointTTF_slow.cc',
  ),
  'orbit_utils': files(
    'core/utils_init.cc',
    'core/type_accessors.cc',
    'utils/wrap_utils.cc',
    'utils/matrix/wrap_matrix.cc',
    'utils/matrix/wrap_phase_vector.cc',
    'utils/wrap_py_base_field_source.cc',
    'utils/wrap_field_source_container.cc',
    'utils/wrap_function.cc',
    'utils/wrap_splinech.cc',
    'utils/statistics/wrap_statmoments2d.cc',
    'utils/bunch/wrap_bunch_extrema_calculator.cc',
    'utils/integration/wrap_gauss_legendre_integrator.cc',
    'utils/polynomial/wrap_polynomial.cc',
    'utils/wrap_numrecipes.cc',
    'utils/bunch/wrap_bunch_utils_functions.cc',
    'utils/harmonic_analysis/wrap_harmonic_data.cc',
    'utils/wrap_random.cc',
  ),
  'aperture': files(
    'core/aperture_init.cc',
    'core/type_accessors.cc',
    'orbit/Apertures/wrap_aperture.cc',
    'orbit/Apertures/wrap_TAperture.cc',
    'orbit/Apertures/wrap_PhaseAperture.cc',
    'orbit/Apertures/wrap_EnergyAperture.cc',
    'orbit/Apertures/wrap_BaseAperture.cc',
    'orbit/Apertures/wrap_PyBaseApertureShape.cc',
    'orbit/Apertures/wrap_PrimitiveApertureShape.cc',
    'orbit/Apertures/wrap_CompositeApertureShape.cc',
    'orbit/Apertures/wrap_ConvexApertureShape.cc',
  ),
  'foil': files(
    'core/foil_init.cc',
    'core/type_accessors.cc',
    'orbit/MaterialInteractions/wrap_foil.cc',
  ),
  'field_sources': files(
    'core/field_sources_init.cc',
    'core/type_accessors.cc',
    'utils/field_sources/wrap_field_sources_module.cc',
    'utils/field_sources/wrap_magnetic_field_source_grid3d.cc',
    'utils/field_sources/wrap_quad_field_source.cc',
    'utils/field_sources/wrap_dipole_field_source.cc',
    'utils/field_sources/wrap_loop_field_source.cc',
  ),
  'rfcavities': files(
    'core/rfcavities_init.cc',
    'core/type_accessors.cc',
    'orbit/RFCavities/wrap_rfcavities.cc',
    'orbit/RFCavities/wrap_Frequency_Cav.cc',
    'orbit/RFCavities/wrap_Harmonic_Cav.cc',
    'orbit/RFCavities/wrap_Barrier_Cav.cc',
    'orbit/RFCavities/wrap_Dual_Harmonic_Cav.cc',
  ),
  'impedances': files(
    'core/impedances_init.cc',
    'core/type_accessors.cc',
    'orbit/Impedances/wrap_impedances.cc',
    'orbit/Impedances/wrap_LImpedance.cc',
    'orbit/Impedances/wrap_TImpedance.cc',
  ),
  'fieldtracker': files(
    'core/fieldtracker_init.cc',
    'core/type_accessors.cc',
    'orbit/FieldTracker/wrap_fieldtracker.cc',
  ),
  'collimator': files(
    'core/collimator_init.cc',
    'core/type_accessors.cc',
    'orbit/MaterialInteractions/wrap_collimator.cc',
  ),
  'error_base': files(
    'core/error_base_init.cc',
    'orbit/Errors/wrap_errorbase.cc',
  ),
}

foreach name, src_list : ext_modules
  python.extension_module(
  name,
  sources: src_list,
  include_directories: inc,
  link_with: core_lib,
  dependencies: dependencies,
  install: true,
  subdir: 'orbit/core',
  )
endforeach
