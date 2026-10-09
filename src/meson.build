# Add Python installation details
python = import('python').find_installation('python3', pure: false)

# Add C++ compiler details
cpp = meson.get_compiler('cpp')

# Add -march=native if the compiler supports it
if cpp.has_argument('-march=native')
  add_project_arguments('-march=native', language: 'cpp')
endif

dependencies = []

dependencies += python.dependency()
dependencies += dependency('fftw3', version: '>= 3.0.0', required: true)

# Detecting if MPICH or OPENMPI are installed and enabling support if present

mpi_use = get_option('USE_MPI')

if mpi_use == 'mpich'
  message('Requested to use MPICH as the MPI implementation.')
  dependencies += dependency('mpich', version: '>= 4.0.0', required: true)
  add_project_arguments('-DUSE_MPI=1', language: 'cpp')

  # Configure dependencies or settings specific to MPICH
elif mpi_use == 'ompi'
  message('Requested to use OpenMPI as the MPI implementation.')
  dependencies += dependency('ompi', version: '>= 4.0.0', required: true)
  add_project_arguments('-DUSE_MPI=1', language: 'cpp')

elif mpi_use == 'none'
  message('Requested to not use MPI.')

else
   mpich_dependency = dependency('mpich', version: '>= 4.0.0', required: false)
   openmpi_dependency = dependency('ompi', version: '>= 4.0.0', required: false)

   if mpich_dependency.found()
     add_project_arguments('-DUSE_MPI=1', language: 'cpp')
     dependencies += mpich_dependency
     message('Using MPICH as the MPI implementation.')
   elif openmpi_dependency.found()
     add_project_arguments('-DUSE_MPI=1', language: 'cpp')
     dependencies += openmpi_dependency
     message('Using OpenMPI as the MPI implementation.')
   else
     message('MPI will not be used.')
   endif

endif

sources = files([
  'linac/wrap_linacmodule.cc',
  'linac/tracking/linac_tracking.cc',
  'linac/tracking/wrap_linac_tracking.cc',
  'linac/rfgap/wrap_BaseRfGap.cc',
  'linac/rfgap/BaseRfGap_slow.cc',
  'linac/rfgap/wrap_RfGapTTF_slow.cc',
  'linac/rfgap/wrap_RfGapThreePointTTF.cc',
  'linac/rfgap/BaseRfGap.cc',
  'linac/rfgap/RfGapThreePointTTF.cc',
  'linac/rfgap/MatrixRfGap.cc',
  'linac/rfgap/wrap_BaseRfGap_slow.cc',
  'linac/rfgap/wrap_RfGapThreePointTTF_slow.cc',
  'linac/rfgap/wrap_SuperFishFieldSource.cc',
  'linac/rfgap/RfGapTTF.cc',
  'linac/rfgap/SuperFishFieldSource.cc',
  'linac/rfgap/RfGapTTF_slow.cc',
  'linac/rfgap/RfGapThreePointTTF_slow.cc',
  'linac/rfgap/wrap_RfGapTTF.cc',
  'linac/rfgap/wrap_MatrixRfGap.cc',
  'orbit/OrbitConst.cc',
  'orbit/wrap_bunch.cc',
  'orbit/wrap_syncpart.cc',
  'orbit/Bunch.cc',
  'orbit/SyncPart.cc',
  'orbit/BunchDiagnostics/wrap_bunch_tune_analysis.cc',
  'orbit/BunchDiagnostics/wrap_bunch_twiss_analysis.cc',
  'orbit/BunchDiagnostics/BunchTwissAnalysis.cc',
  'orbit/BunchDiagnostics/BunchTuneAnalysis.cc',
  'orbit/FieldTracker/wrap_fieldtracker.cc',
  'orbit/FieldTracker/FieldTracker.cc',
  'orbit/MaterialInteractions/wrap_foil.cc',
  'orbit/MaterialInteractions/Foil.cc',
  'orbit/MaterialInteractions/Collimator.cc',
  'orbit/MaterialInteractions/numrecipes.cc',
  'orbit/MaterialInteractions/cross_sections.cc',
  'orbit/MaterialInteractions/wrap_collimator.cc',
  'orbit/MaterialInteractions/MaterialInteractions.cc',
  'orbit/RFCavities/Dual_Harmonic_Cav.cc',
  'orbit/RFCavities/Frequency_Cav.cc',
  'orbit/RFCavities/Barrier_Cav.cc',
  'orbit/RFCavities/wrap_Harmonic_Cav.cc',
  'orbit/RFCavities/Harmonic_Cav.cc',
  'orbit/RFCavities/wrap_Dual_Harmonic_Cav.cc',
  'orbit/RFCavities/wrap_Barrier_Cav.cc',
  'orbit/RFCavities/wrap_Frequency_Cav.cc',
  'orbit/RFCavities/wrap_rfcavities.cc',
  'orbit/Impedances/wrap_TImpedance.cc',
  'orbit/Impedances/TImpedance.cc',
  'orbit/Impedances/wrap_impedances.cc',
  'orbit/Impedances/LImpedance.cc',
  'orbit/Impedances/wrap_LImpedance.cc',
  'orbit/SynchPartRedefinition/SynchPartRedefinitionZdE.cc',
  'orbit/SynchPartRedefinition/wrap_synch_part_redefinition_z_de.cc',
  'orbit/Errors/errorbase.cc',
  'orbit/Errors/wrap_errorbase.cc',
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
  'orbit/Apertures/wrap_aperture.cc',
  'orbit/Apertures/wrap_CompositeApertureShape.cc',
  'orbit/Apertures/RectangularApertureShape.cc',
  'orbit/Apertures/PyBaseApertureShape.cc',
  'orbit/Apertures/ConvexApertureShape.cc',
  'orbit/Apertures/wrap_ConvexApertureShape.cc',
  'orbit/Apertures/wrap_PyBaseApertureShape.cc',
  'orbit/Apertures/EnergyAperture.cc',
  'orbit/Apertures/BaseAperture.cc',
  'orbit/Apertures/EllipseApertureShape.cc',
  'orbit/Apertures/wrap_BaseAperture.cc',
  'orbit/Apertures/wrap_PrimitiveApertureShape.cc',
  'orbit/Apertures/CompositeApertureShape.cc',
  'orbit/Apertures/wrap_EnergyAperture.cc',
  'orbit/Apertures/wrap_TAperture.cc',
  'orbit/Apertures/wrap_PhaseAperture.cc',
  'orbit/Apertures/BaseApertureShape.cc',
  'trackerrk4/wrap_runge_kutta_tracker.cc',
  'trackerrk4/wrap_py_external_effects.cc',
  'trackerrk4/ExternalEffects.cc',
  'trackerrk4/ExtEffectsContainer.cc',
  'trackerrk4/wrap_trackerrk4.cc',
  'trackerrk4/PyExternalEffects.cc',
  'trackerrk4/RungeKuttaTracker.cc',
  'trackerrk4/wrap_ext_effects_container.cc',
  'mpi/wrap_orbit_mpi.cc',
  'mpi/wrap_mpi_comm.cc',
  'mpi/orbit_mpi.cc',
  'mpi/wrap_mpi_datatype.cc',
  'mpi/wrap_mpi_request.cc',
  'mpi/wrap_mpi_op.cc',
  'mpi/wrap_mpi_status.cc',
  'mpi/wrap_mpi_group.cc',
  'spacecharge/PoissonSolverFFT3D.cc',
  'spacecharge/SpaceChargeCalcUnifEllipse.cc',
  'spacecharge/ForceSolver2D.cc',
  'spacecharge/PoissonSolverFFT2D.cc',
  'spacecharge/wrap_poissonsolverfft2d.cc',
  'spacecharge/wrap_poissonsolverfft3d.cc',
  'spacecharge/PoissonSolver2D.cc',
  'spacecharge/SpaceChargeForceCalc2p5D.cc',
  'spacecharge/PoissonSolver3D.cc',
  'spacecharge/SpaceChargeCalc3D.cc',
  'spacecharge/wrap_boundary2d.cc',
  'spacecharge/wrap_forcesolverfft2d.cc',
  'spacecharge/ForceSolverFFT2D.cc',
  'spacecharge/SpaceChargeCalcSliceBySlice2D.cc',
  'spacecharge/Grid3D.cc',
  'spacecharge/wrap_spacecharge.cc',
  'spacecharge/wrap_spacechargecalc2p5d_rb.cc',
  'spacecharge/Grid2D.cc',
  'spacecharge/SpaceChargeCalc2p5Drb.cc',
  'spacecharge/Grid1D.cc',
  'spacecharge/LSpaceChargeCalc.cc',
  'spacecharge/wrap_grid1D.cc',
  'spacecharge/wrap_uniform_ellipsoid_field_calculator.cc',
  'spacecharge/ShapedBoundary2D.cc',
  'spacecharge/wrap_grid3D.cc',
  'spacecharge/BaseBoundary2D.cc',
  'spacecharge/wrap_grid2D.cc',
  'spacecharge/wrap_spacechargecalc_uniform_ellipse.cc',
  'spacecharge/wrap_spacechargecalc3d.cc',
  'spacecharge/UniformEllipsoidFieldCalculator.cc',
  'spacecharge/wrap_spacechargeforcecalc2p5d.cc',
  'spacecharge/SpaceChargeCalc2p5D.cc',
  'spacecharge/wrap_spacechargecalc_slicebyslice_2D.cc',
  'spacecharge/wrap_spacechargecalc2p5d.cc',
  'spacecharge/wrap_lspacechargecalc.cc',
  'utils/wrap_function.cc',
  'utils/FieldSourceContainer.cc',
  'utils/AttributesBucket.cc',
  'utils/bessel.cc',
  'utils/PyBaseFieldSource.cc',
  'utils/wrap_field_source_container.cc',
  'utils/OU_SplineCH.cc',
  'utils/CppPyWrapper.cc',
  'utils/wrap_random.cc',
  'utils/wrap_utils.cc',
  'utils/StringUtils.cc',
  'utils/wrap_py_base_field_source.cc',
  'utils/wrap_splinech.cc',
  'utils/Random.cc',
  'utils/wrap_numrecipes.cc',
  'utils/BaseFieldSource.cc',
  'utils/OU_Function.cc',
  'utils/MathPolynomial.cc',
  'utils/BufferStore.cc',
  'utils/statistics/wrap_statmoments2d.cc',
  'utils/statistics/StatMoments2D.cc',
  'utils/field_sources/wrap_loop_field_source.cc',
  'utils/field_sources/LoopFieldSource.cc',
  'utils/field_sources/QuadFieldSource.cc',
  'utils/field_sources/ShiftedFieldSource.cc',
  'utils/field_sources/wrap_dipole_field_source.cc',
  'utils/field_sources/wrap_quad_field_source.cc',
  'utils/field_sources/wrap_magnetic_field_source_grid3d.cc',
  'utils/field_sources/DipoleFieldSource.cc',
  'utils/field_sources/wrap_field_sources_module.cc',
  'utils/field_sources/MagnetFieldSourceGrid3D.cc',
  'utils/integration/gauss_legendre_points.cc',
  'utils/integration/wrap_gauss_legendre_integrator.cc',
  'utils/integration/GaussLegendreIntegrator.cc',
  'utils/bunch/ParticlesWithIdFunctions.cc',
  'utils/bunch/wrap_bunch_extrema_calculator.cc',
  'utils/bunch/TwissFilteringFunctions.cc',
  'utils/bunch/InitialCoordsAttrFunctions.cc',
  'utils/bunch/BunchExtremaCalculator.cc',
  'utils/bunch/wrap_bunch_utils_functions.cc',
  'utils/harmonic_analysis/HarmonicData.cc',
  'utils/harmonic_analysis/wrap_harmonic_data.cc',
  'utils/polynomial/wrap_polynomial.cc',
  'utils/polynomial/OU_Polynomial.cc',
  'utils/matrix/wrap_matrix.cc',
  'utils/matrix/PhaseVector.cc',
  'utils/matrix/MatrixOperations.cc',
  'utils/matrix/wrap_phase_vector.cc',
  'utils/matrix/Matrix.cc',
  'teapot/wrap_teapotbase.cc',
  'teapot/wrap_matrix_generator.cc',
  'teapot/teapotbase.cc',
  'teapot/MatrixGenerator.cc'
])

inc = include_directories([
  'main',
  'orbit/SynchPartRedefinition',
  'trackerrk4',
  'utils/ellipticalint',
  'teapot',
  'mpi',
  'orbit/ParticlesAttributes',
  'orbit/Errors',
  'orbit/FieldTracker',
  'utils/statistics',
  'utils',
  'linac/tracking',
  'spacecharge',
  'orbit/MaterialInteractions',
  'utils/field_sources',
  'utils/bunch',
  'orbit/Impedances',
  'utils/harmonic_analysis',
  'linac',
  'utils/polynomial',
  'linac/rfgap',
  'orbit/RFCavities',
  'utils/matrix',
  'orbit/BunchDiagnostics',
  'orbit',
  'utils/integration',
  'orbit/Apertures'
])

core_lib = library(
  'core',
  sources: sources,
  include_directories: inc,
  override_options: ['b_lundef=false'],
  dependencies: dependencies,
  install: true,
)

core_dep = declare_dependency(link_with : core_lib)

ext_modules = {
  'orbit_mpi':     'mpi_init.cc',
  'bunch':         'bunch_init.cc',
  'spacecharge':   'spacecharge_init.cc',
  'trackerrk4':    'trackerrk4_init.cc',
  'teapot_base':   'teapot_base_init.cc',
  'linac':         'linac_init.cc',
  'orbit_utils':   'utils_init.cc',
  'aperture':      'aperture_init.cc',
  'foil':          'foil_init.cc',
  'field_sources': 'field_sources_init.cc',
  'rfcavities':    'rfcavities_init.cc',
  'impedances':    'impedances_init.cc',
  'fieldtracker':  'fieldtracker_init.cc',
  'collimator':    'collimator_init.cc',
  'error_base':    'error_base_init.cc',
}

base = meson.current_source_dir() + '/core'

foreach name, src : ext_modules
  python.extension_module(
  name,
  sources: base / src,
  include_directories: inc,
  dependencies: [core_dep],
  install: true,
  subdir: 'orbit/core',
  )
endforeach
