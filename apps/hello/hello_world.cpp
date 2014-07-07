#include <graphlab.hpp>

int main(int argc, char** argv)
{
  //Graphlab initialization
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;

  //Main body
  dc.cout() << "Hello World!\n";


  //Main body ends
  graphlab::mpi_tools::finalize();
  return 0;
}
