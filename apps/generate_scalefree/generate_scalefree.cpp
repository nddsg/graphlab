#include <graphlab.hpp>

int main(int argc, char** argv) {
  // Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;
  global_logger().set_log_level(LOG_INFO);

  size_t vertices = 0;
  std::string outgraph, outformat;

  bool gzip = false;
  bool in_degree = true;
  double alpha = 2.1;
  // Parse command line options -----------------------------------------------
  graphlab::command_line_options clopts("Graph Format Conversion.", true);

  clopts.attach_option("vertices", vertices,
                       "Generates a synthetic powerlaw graph with this many vertices.");
  clopts.attach_option("outgraph", outgraph,
                       "The output graph file. Required ");
  clopts.attach_option("outformat", outformat,
                       "The output graph file format");
  clopts.attach_option("outgzip", gzip,
                       "If output is to be gzip compressed");
  clopts.attach_option("in_degree", in_degree,
                       "If true, the graph will have power-law in-degree, else uniform.");
  clopts.attach_option("alpha", alpha,
                       "Alpha parameter in the power law distribution.");

  if(!clopts.parse(argc, argv)) {
    dc.cout() << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }
  if (vertices==0 || outgraph.length() == 0) {
    clopts.print_description();
    return EXIT_FAILURE;
  }
  typedef graphlab::distributed_graph<graphlab::empty, graphlab::empty> graph_type;
  graph_type graph(dc, clopts);

  graph.load_synthetic_powerlaw(vertices, in_degree, alpha, 100000000 /*max degree*/);

  graph.finalize();

  dc.cout() << "#vertices: " << graph.num_vertices()
            << " #edges:" << graph.num_edges() << std::endl;

  graph.save_format(outgraph, outformat, gzip);

  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // End of main
