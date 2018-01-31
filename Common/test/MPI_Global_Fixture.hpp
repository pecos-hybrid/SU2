#include <boost/test/included/unit_test.hpp>

#include "../../Common/include/config_structure.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

struct MPIGlobalFixture {

  MPIGlobalFixture() {
#ifdef HAVE_MPI
  MPI_Init(NULL,NULL);
  BOOST_TEST_MESSAGE("Calling MPI_Init...");
#endif
  }

  ~MPIGlobalFixture() {
#ifdef HAVE_MPI
  MPI_Finalize();
  BOOST_TEST_MESSAGE("Calling MPI_Finalize...");
#endif
  }

};

BOOST_GLOBAL_FIXTURE( MPIGlobalFixture );
