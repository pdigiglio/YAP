set(INCLUDE_DIRECTORIES
	${YAP_SOURCE_DIR}/include)
include_directories(${INCLUDE_DIRECTORIES})

set(YAP_SOURCES
	BlattWeisskopf.cxx
	BreitWigner.cxx
	CachedDataValue.cxx
	CachedValue.cxx
	ClebschGordan.cxx
	DataAccessor.cxx
	DataPartition.cxx
	DataPoint.cxx
	DataSet.cxx
	DecayChannel.cxx
	DecayTree.cxx
	DecayingParticle.cxx
	FinalStateParticle.cxx
	Flatte.cxx
	FreeAmplitude.cxx
	FourMomenta.cxx
	HelicityAngles.cxx
	HelicityFormalism.cxx
  ImportanceSampler.cxx
	MassShape.cxx
	MassShapeWithNominalMass.cxx
	MeasuredBreakupMomenta.cxx
	Model.cxx
  ModelIntegral.cxx
	ParticleCombination.cxx
	ParticleCombinationCache.cxx
	Particle.cxx
	ParticleFactory.cxx
	PoleMass.cxx
	QuantumNumbers.cxx
	RecalculableDataAccessor.cxx
	Resonance.cxx
	SpinAmplitude.cxx
	StatusManager.cxx
	WignerD.cxx
	ZemachFormalism.cxx
)

add_library(YAP SHARED ${YAP_SOURCES})

# install destinations can be passed via the command line:
# cmake -DLIBRARY_OUTPUT_DIRECTORY:PATH=<lib_path> <path_to_CMakeLists.tex>
# otherwise, default LD_LIBRARY_PATH
if(NOT DEFINED LIBRARY_OUTPUT_DIRECTORY)
	set(LIBRARY_OUTPUT_DIRECTORY ${YAP_SOURCE_DIR}/lib/${CMAKE_BUILD_TYPE})
endif()

if(NOT DEFINED INCLUDE_OUTPUT_DIRECTORY)
	set(INCLUDE_OUTPUT_DIRECTORY ${YAP_SOURCE_DIR}/include/YAP)
endif()

install(TARGETS YAP LIBRARY DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})

# Matches all the headers in ${YAPDID}/include and its subdirs
file(GLOB_RECURSE
	 INSTALL_INCLUDES ${YAP_SOURCE_DIR}/include/*.h)

#message(STATUS "${INSTALL_INCLUDES}")

#install(FILE       ${INSTALL_INCLUDES}
#	    DESTINATION ${INCLUDE_OUTPUT_DIRECTORY}
#	)
