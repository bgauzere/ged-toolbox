#include <vector>

#ifndef __PGRAPHH__
#define __PGRAPHH__

class GEdge{
private:
   /** The next neighbourhood node. */
   GEdge *next; 
   /** The rank of the node in the graph array. */
   long incident_node; //GNode ?
   /** The number which identifies the object. */
   long edge_id;

public:
   /** The attribute of the edge */
   int attr;
  
   /**
    * Creates a new edge to the specified node, with
    * the specified weight and the specified next edge.
    * @param n	the connected node.
    * @param adj the next edge.
    * @param attr the label.
    */
   GEdge( long n, GEdge *adj, int attr ): next(adj), incident_node(n), edge_id(-1), attr(attr) {};

   /**
    * Creates a new edge to the specified incident_node, with
    * the specified weight and the specified next edge.
    * @param n	the connected incident_node.
    * @param adj	the next edge.
    * @param i       the identifier of the edge
    * @param w	the weight.
    */
   GEdge( long n, GEdge *adj, long i,int attr): next(adj), incident_node(n), edge_id(i), attr(attr) {};

   /**
    * Deletes the edge
    */
   ~GEdge(){};

   /**
    * Returns the number of the incident incident_node.
    * @return	the number of the connected incident_node.
    */
   long IncidentNode() const { return incident_node; };

   /**
    * Returns the next neighbourhood incident_node.
    * @return	the next edge.
    */
   GEdge* Next() const { return next; };

   /**
    * Sets a new neighbourhood incident_node.
    * @param n	the  neighbourhood edge.
    * @return	the neighbourhood edge.
    */
   GEdge* Next( GEdge* n ) { return next=n; };

   /**
    * Returns the index of the referenced object.
    * @return	the index.
    */
   long EdgeId() const { return edge_id; }

   /**
    * Sets the new index of the referenced object.
    * @param i	the new index.
    * @return	the new index of the object.
    */
   long EdgeId( long i ) { return edge_id=i; }

};

/** @brief A node of a graph.
 *
 * The class <code>GNode</code> defines a node.
 * A node indexes an object in a separate array of objects
 * at specified coordinates in the related image.
 * It is characterized by a value and a list
 * of connected nodes.
 */
class GNode{
private:
  /** The list of the adjacent nodes. */
  GEdge *adjacents;
  /** The number which identifies the object. */
  long item;
  // /** Trash for removed edges */
  // std::vector<GEdge *> etrash;  //XXX : A virer ?
  
public :
  /** The valuation of the node. */
  int attr;	

  /**
   * Creates a new node with the specified index,
   * and the specified coordinates.
   * @param i	the index of the referenced object.
   * @param p	the specified coordinates.
   */
  GNode( long i, int attr ): adjacents(0), attr(attr) { };

  /**
   * Deletes the node.
   */
  ~GNode();

  /**
   * Returns the list of all the connected nodes.
   * @return	the list of connected nodes.
   */
  GEdge* getIncidentEdges() const { return adjacents; };

  GEdge * Connect( long incidentNode, int attr);
  GEdge * Connect( long incidentNode, long edge_id, int attr);

  int Degree();
  /**
   * Deletes the specified node from the list of connected nodes.
   * @param n	the specified node.
   * @return the new list of edges.
   */
  GEdge* UnConnect( long n );

  /**
   * Returns the index of the referenced object.
   * @return	the index.
   */
  long Item() const { return item; }

  /**
   * Sets the new index of the referenced object.
   * @param i	the new index.
   * @return	the new index of the object.
   */
  long Item( long i ) { return item=i; }
};

/** @brief A 2D graph.
 *
 * A graph is a set of nodes connected to some other nodes.
 * The two types of graph are supported: directed and undirected;
 * the type must be specified with the constructor.
 * A node is characterized by a value which determines if the
 * node is active in the representation of the graph or not.
 * A node can be any of the Pandore object. It indexes an 
 * item in a separate array of items. (The secret is to use an
 * item as a pointer to a specific objet in an array [no type]).
 * <br>For the use of Graph2d see @ref graph_page.
 */
class Graph {
private :
  std::vector<GNode *> tnode;
  long nbNodes;
  long nbEdges;
  bool _directed;

  friend class GEdge;
   
public :
  /**
   * Deletes the graph.
   */
  ~Graph();

  /**
   * @return true if the graph is directed.
   */
  bool isDirected() const { return _directed; }
  /**
   * Returns the number of nodes.
   * @return	the size.
   */
  long Size() const { return nbNodes; }
  long getNbEdges() const { return _directed?nbEdges:nbEdges/2; }
  /**
   * Creates a new graph with no data.
   * @param directed true for creating a directed graph.
   */
  Graph( bool directed =false): tnode(0), nbNodes(0), nbEdges(0), _directed(directed) { }
  /**
   * Initialize a graph nb_nodes and corresponding to adjacency matrix am
   * the format of am is the following : 
   * am[i][j] : label of edge between nodes i and j, 0 if no edge
   * am[i][i] : label of node i
   */
  Graph(int * am, int nb_nodes, bool directed=false);

  /**
   * Returns the node at the specified coordinates.
   * @param pos	the coordinates.
   * @return	the node at the specified coordinates.
   */
  GNode *operator[]( long pos ){ return(tnode[pos]); }

  /**
   * Returns the node at the specified coordinates.
   * @param pos	the coordinates.
   * @return	the node at the specified coordinates.
   */
  const GNode *operator[]( long pos ) const { return(tnode[pos]); }
   
  /**
   * Adds the specified node that references the specified item. 
   * @param node	the new Node
   * @return	index du noeud.
   */
  int Add( GNode* node );
      
  /**
   * Deletes the specified node from the graph. Unlinks it
   * from connected nodes.
   * @param s	the node to be deleted.
   * @return	SUCCESS or FAILURE.
   */
  GNode * Del( long s );

  GEdge * Link(long firstNode, long SecondNode , int label);
  //TODO : UnLink !
  bool isLinked(long firstNode, long SecondNode);
  GEdge * getEdge(long firstNode, long SecondNode);
  GEdge * getSymmetricEdge(long nodeId,const GEdge * p);
};

#endif // __PGRAPHH__
