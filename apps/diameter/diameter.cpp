#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <boost/unordered_map.hpp>
#include <set>
#include <iostream>
#include <ctime>
#include <graphlab/serialization/oarchive.hpp>
#include "BitPackMap32.hpp"
#include "VertexHistoryPack.hpp"

// We must include this std::map.operator+= definition before including graphlab.hpp
template <typename Key, typename T, typename Comp, typename Allocator>
inline std::map< Key, T, Comp, Allocator >& operator+=(std::map< Key, T, Comp, Allocator >& lvalue, const std::map< Key, T, Comp, Allocator >& rvalue) {
  for (typename std::map< Key, T, Comp, Allocator >::const_iterator it = rvalue.begin(); it != rvalue.end(); ++it) {
    lvalue[it->first] += it->second;
  }
  return lvalue;
} // end of operator +=

template <typename Key, typename T, typename Comp, typename Allocator>
inline boost::unordered_map< Key, T, Comp, Allocator >& operator+=(boost::unordered_map< Key, T, Comp, Allocator >& lvalue, const boost::unordered_map< Key, T, Comp, Allocator >& rvalue) {
  if (lvalue.size() + rvalue.size() > lvalue.load_factor() * lvalue.bucket_count()) {
    lvalue.reserve(lvalue.size() + rvalue.size());
  }
  for (typename boost::unordered_map< Key, T, Comp, Allocator >::const_iterator it = rvalue.begin(); it != rvalue.end(); ++it) {
    lvalue[it->first] += it->second;
  }
  return lvalue;
} // end of operator +=

// std::vector.operator+= definition
template <typename T>
inline std::vector< T >& operator+=(std::vector< T >& lvalue, const std::vector< T >& rvalue) {
  if(!rvalue.empty()) {
    if(lvalue.empty()) lvalue = rvalue;
    else {
      for(size_t t = 0; t < lvalue.size(); ++t) lvalue[t] += rvalue[t];
    }
  }
  return lvalue;
} // end of operator +=

// std::set.operator+= definition
template <typename T, typename Comp, typename Allocator>
inline std::set< T, Comp, Allocator >& operator+=(std::set< T, Comp, Allocator >& lvalue, const std::set< T, Comp, Allocator >& rvalue) {
  lvalue.insert(rvalue.begin(), rvalue.end());
  return lvalue;
} // end of operator +=

#include <graphlab.hpp>

// Globals
bool IS_DIRECTED = false;
uint32_t DIA_ARR_SIZE = 0;
uint32_t BATCHES = 1;

// Type Definitions
typedef BitPackMap vertex_map;
//typedef BitPackMap vertex_vector;
typedef std::map<uint32_t, uint32_t> uint32_map;


struct VertexHistory {
  BitPackMap m;
  // std::vector<uint32_t> v;
  VertexHistory() {
    clear();
  }
  inline VertexHistory& operator+=(VertexHistory& rvalue) {
    m += rvalue.m;
    // v += rvalue.v;
    return *this;
  }
  // uint32_t & operator[] (const graphlab::vertex_id_type & key) {
    // return v[key];
  // }
  bool get(uint32_t pos) {
    return m[pos];
  }
  void set(uint32_t pos, bool value) {
    m.set(pos, value);
  }
  void clear() {
    m.clear();
    // v.reserve(DIA_ARR_SIZE);
    // for (unsigned int i = 0; i < DIA_ARR_SIZE; i++) {
      // v[i] = 0;
    // }
  }
  void save(graphlab::oarchive& oarc) const {
    // oarc << m << v;
  }
  void load(graphlab::iarchive& iarc) {
    // iarc >> m >> v;
  }
};


graphlab::vertex_id_type id2offset(graphlab::vertex_id_type id) {
  return id / BATCHES;
}

template <typename T, typename Comp, typename Allocator>
inline bool in_set(const std::set< T, Comp, Allocator > & s, const uint32_t val) {
  return s.find(val) != s.end();
}

struct vdata {
  //vertex_map distance;
  //uint32_map distance_count;
  bool active;
  vertex_map juggle[4];
  VertexHistoryPack history;
  int sent_current;
  int sent_previous;
  //int sent_old;
  int temp;
  int max_distance;
  // Constructor
  vdata() : active(false), sent_current(0), sent_previous(1), /*sent_old(2), */temp(3), max_distance(0), history() {
  }
  // Copy Constructor
/*  vdata(const vdata & other) : history() {
    //distance = other.distance;
    //distance_count = other.distance_count;
    active = other.active;
    sent_current = other.sent_current;
    sent_previous = other.sent_previous;
    sent_old = other.sent_old;
    temp = other.temp;
    history = other.history;
    max_distance = other.max_distance;
  }
*/  // Deconstructor
  ~vdata() {
  }
  // Assignment Operator
/*  vdata& operator=(const vdata& rhs) {
    //distance = rhs.distance;
    //distance_count = rhs.distance_count;
    active = rhs.active;
    sent_current = rhs.sent_current;
    sent_previous = rhs.sent_previous;
    sent_old = rhs.sent_old;
    temp = rhs.temp;
    max_distance = rhs.max_distance;
    if (IS_DIRECTED) {
      //history.reserve(DIA_ARR_SIZE);
      history.clear();
    }
    history = rhs.history;
    return *this;
  }
*/  void clear_juggle2() {
    history.clear();
  }
  void save(graphlab::oarchive& oarc) const {
    oarc << sent_current << sent_previous << /*distance << distance_count <<*/ active << /*sent_old << *//*juggle[0] << juggle[1] << juggle[2] << juggle[3] << */history;
  }
  void load(graphlab::iarchive& iarc) {
    iarc >> sent_current >> sent_previous >> /*distance >> distance_count >>*/ active >> /*sent_old >> *//*juggle[0] >> juggle[1] >> juggle[2] >> juggle[3] >> */history;
  }
};

typedef graphlab::distributed_graph<vdata, graphlab::empty> graph_type;

class diameter: 
  public graphlab::ivertex_program<graph_type, vertex_map> {
public:
  edge_dir_type gather_edges(icontext_type& context, const vertex_type& vertex) const {
    return IS_DIRECTED ? graphlab::IN_EDGES : graphlab::ALL_EDGES;
  }

  vertex_map gather(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
    return (vertex.id() == edge.source().id()) ? edge.target().data().juggle[edge.target().data().sent_current] : edge.source().data().juggle[edge.source().data().sent_current];
  }

  void apply(icontext_type& context, vertex_type& vertex, const gather_type& incoming) {
    vertex_data_type & data = vertex.data();

    if (data.active && incoming.empty()) {
      // This is the initialization message.
      // Setup the vertex's messaging.
      uint32_t id = vertex.id();
      data.juggle[data.temp].insert(id2offset(id));
      data.active = false;
      // std::cout << "Activating " << vertex.id() << " with message #" << id2offset(id) << "\n";
    }
    
    //uint32_t count = 0;
    uint32_t new_temp;

    // Directed Graph
    if (IS_DIRECTED) {
      //data.juggle[data.temp].reserve(incoming.size() + 1);
      // Move the incoming into the temp map
      for (vertex_map::iterator it = incoming.begin(); it != incoming.end(); ++it) {
        // std::cout << "\t" << vertex.id() << " receiving message from " << it.first << "\n";
        if (!data.history.get(it.first)) {
          data.juggle[data.temp].insert(it.first);
          data.max_distance = std::max(data.max_distance, context.iteration());
          // std::cout << "\t\tPassing along\n";
        }
        else {
          // std::cout << "\t\tIgnoring\n";
        }
        //if ((it->first != vertex.id()) /*&& (data.distance[it->first] < 1)*/) {
          //data.distance[it->first] = context.iteration();
          //count++;
        //}
      }
      for (vertex_map::iterator it = data.juggle[data.temp].begin(); it != data.juggle[data.temp].end(); ++it) {
        //data.history.insert(it.first);
        data.history.set(it.first, true);
      }
      // In a directed graph, temp is moved to sent_current, and the old
      // sent_current is recycled for the next iteration's temp.
      new_temp = data.sent_current;
      data.sent_current = data.temp;
      data.temp = new_temp;
      data.juggle[new_temp].clear();
    }
    // Undirected Graph
    else {
      for (vertex_map::iterator it = incoming.begin(); it != incoming.end(); ++it) {
        if (!data.juggle[data.sent_previous][it.first]
          && !data.juggle[data.sent_current][it.first]
          /*&& !data.juggle[data.sent_old][it.first]*/) {
          data.juggle[data.temp].insert(it.first);
          data.max_distance = std::max(data.max_distance, context.iteration());
        }
        //if ((it->first != vertex.id()) /*&& (data.distance[it->first] < 1)*/) {
          //data.distance[it->first] = context.iteration();
          //count++;
        //}
      }

      // In an undirected graph, sent_old is discarded during the variable
      // rotation, so it will become the new temp variable.
      new_temp = /*data.sent_old;
      data.sent_old = */data.sent_previous;
      data.sent_previous = data.sent_current;
      data.sent_current = data.temp;
      data.temp = new_temp;
      data.juggle[new_temp].clear();
    }

    //data.distance_count[context.iteration()] = count;
  }

  edge_dir_type scatter_edges(icontext_type& context, const vertex_type& vertex) const {
    if (!vertex.data().juggle[vertex.data().sent_current].empty()) {
      return IS_DIRECTED ? graphlab::OUT_EDGES : graphlab::ALL_EDGES;
    }
    return graphlab::NO_EDGES;
  }

  void scatter(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
    context.signal(vertex.id() != edge.target().id() ? edge.target() : edge.source());
  }

  void save(graphlab::oarchive& oarc) const {
    //oarc << last_change;
  }

  void load(graphlab::iarchive& iarc) {
    //iarc >> last_change;
  }
};

struct MaxID {
  graphlab::vertex_id_type id;
  MaxID() : id(0) {}
  MaxID(graphlab::vertex_id_type m) : id(m) {}
  // operator +=
  MaxID & operator+= (const MaxID &rhs) {
    id = std::max(id, rhs.id);
    return *this;
  }
  void save(graphlab::oarchive& oarc) const {
    oarc << id;
  }
  void load(graphlab::iarchive& iarc) {
    iarc >> id;
  }
};

MaxID find_maxid(graphlab::synchronous_engine<diameter>::icontext_type& context, const graph_type::vertex_type& vertex) {
  return MaxID(vertex.id());
}

struct Maxpath {
  uint32_t length;
  Maxpath() : length(0) {}
  Maxpath(uint32_t m) : length(m) {}
  // operator +=
  Maxpath & operator+= (const Maxpath &rhs) {
    length = std::max(length, rhs.length);
    return *this;
  }
  void save(graphlab::oarchive& oarc) const {
    oarc << length;
  }
  void load(graphlab::iarchive& iarc) {
    iarc >> length;
  }
};

Maxpath find_diameter(graphlab::synchronous_engine<diameter>::icontext_type& context, const graph_type::vertex_type& vertex) {
  //uint32_t length = 0;
  return Maxpath(vertex.data().max_distance);/*
  for (vertex_map::const_iterator it = vertex.data().distance.begin(); it != vertex.data().distance.end(); ++it) {
    length = std::max(length, it->second);
  }
  return Maxpath(length);*/
}

class Distance_dump {
  uint32_t diameter;
  public:
  Distance_dump(uint32_t d) : diameter(d) {}
  std::string save_vertex(graph_type::vertex_type v) {
    std::stringstream strm;
    strm << v.id();
    //for (uint32_t i = 0; i < diameter; i++) {
    //  strm << "\t" << v.data().distance_count[i + 1];
    //}
    strm << "\n";
    return strm.str();
  }
  std::string save_edge(graph_type::edge_type e) { return ""; }
};

class Activate_vertex {
  static uint32_t current, total;
  public:
  Activate_vertex(uint32_t c, uint32_t t) {
    current = c;
    total = t;
  }
  static void activate(graphlab::synchronous_engine<diameter>::icontext_type& context, graph_type::vertex_type vertex) {
    vertex.data().active = (vertex.id() % total == current ? true : false);
    if (IS_DIRECTED) {
      vertex.data().clear_juggle2();
      vertex.data().history.resize(DIA_ARR_SIZE);
    }
    else {
      vertex.data().juggle[vertex.data().sent_previous].clear();
      //vertex.data().juggle[vertex.data().sent_old].clear();
    }
    vertex.data().juggle[vertex.data().sent_current].clear();
    vertex.data().juggle[vertex.data().temp].clear();
  }
};
uint32_t Activate_vertex::current = 0, Activate_vertex::total = 0;

std::string getTime(clock_t ticks) {
  int permin = 60;
  int perhour = 60 * permin;
  float seconds = ((float)ticks) / (float)1000;
  int hours = seconds / perhour;
  seconds -= hours * perhour;
  int minutes = seconds / permin;
  seconds -= minutes * permin;
  char buff[100];
  sprintf(buff, "%d:%02d:%05.2f", hours, minutes, seconds);
  return std::string(buff);
}

int main(int argc, char** argv) {
  std::cout << "Graph diameter\n\n";
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;

  std::string datafile;
  //parse command line
  graphlab::command_line_options clopts("Graph diameter.");
  std::string graph_dir;
  std::string format = "adj";
  clopts.attach_option("graph", graph_dir,
                       "The graph file. This is not optional");
  clopts.add_positional("graph");
  clopts.attach_option("format", format,
                       "The graph file format");
  clopts.attach_option("batches", BATCHES, "Process the graph in successive groups, default 1.");

  clopts.attach_option("directed", IS_DIRECTED, "Whether or not this is a directed graph.");

  if (!clopts.parse(argc, argv)){
    dc.cout() << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }
  if (graph_dir == "") {
    std::cout << "--graph is not optional\n";
    return EXIT_FAILURE;
  }

  // Load graph
  graph_type graph(dc, clopts);
  dc.cout() << "Loading graph in format: " << format << std::endl;
  graph.load_format(graph_dir, format);
  graph.finalize();

  graphlab::timer timer;
  double itotal, iend;
  clock_t time;
  time = clock();
  
  // Initialize vertices
  graphlab::synchronous_engine<diameter> engine(dc, graph, clopts);

  // Main iteration
  MaxID maxid = engine.map_reduce_vertices<MaxID>(find_maxid);
  dc.cout() << "Max ID: " << maxid.id << "\n";
  DIA_ARR_SIZE = ceil((maxid.id + 1) / (float)BATCHES);
  dc.cout() << "Array size: " << DIA_ARR_SIZE << "\n";
  
  itotal = 0;
  for (uint32_t i = 0; i < BATCHES; i++) {
    dc.cout() << "Beginning batch " << (i + 1) << " of " << BATCHES << "\n";
    timer.start();
    engine.transform_vertices(Activate_vertex(i, BATCHES).activate);
    engine.signal_all();
    engine.start();
    iend = timer.current_time_millis();
    itotal += (iend);
    double projected = itotal * BATCHES / (float) (i + 1);
    dc.cout() << "\tIteration time:\t" << getTime(iend) << "\n\tTotal spent:\t" << getTime(itotal) << "\n\tProjected Total Time:\t\t" << getTime(projected) << "\n\tProjected Time Remaining:\t" << getTime(projected - itotal) << "\n"; 
  }
  Maxpath m = engine.map_reduce_vertices<Maxpath>(find_diameter);

  dc.cout() << "TIME: " << (itotal / (float) 1000) << "\n";
  dc.cout() << "CLOCKTIME: " << ((float)(clock() - time) / CLOCKS_PER_SEC) << "\n";
  dc.cout() << "Diameter: " << m.length << "\n";

/*
  graph.save("distance_count",
    Distance_dump(m.length),
    false, // set to true if each output file is to be gzipped
    true, // whether vertices are saved
    false); // whether edges are saved
*/
  graphlab::mpi_tools::finalize();

  return EXIT_SUCCESS;
}
