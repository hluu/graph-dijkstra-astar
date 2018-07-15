/**
 * @author UCSD MOOC development team and YOU
 * 
 * A class which reprsents a graph of geographic locations
 * Nodes in the graph are intersections between 
 *
 */
package roadgraph;

import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;
import java.util.Stack;
import java.util.function.Consumer;

import geography.GeographicPoint;
import util.GraphLoader;

/**
 * @author UCSD MOOC development team and YOU
 * 
 *         A class which represents a graph of geographic locations Nodes in the
 *         graph are intersections between
 *
 */
public class MapGraph {
	// TODO: Add your member variables here in WEEK 3
	private Map<GeographicPoint, List<StreetSegment>> graph;

	private int numEdges;

	/**
	 * Create a new empty MapGraph
	 */
	public MapGraph() {
		// TODO: Implement in this constructor in WEEK 3
		graph = new HashMap<>();
	}

	/**
	 * Get the number of vertices (road intersections) in the graph
	 * 
	 * @return The number of vertices in the graph.
	 */
	public int getNumVertices() {
		// TODO: Implement this method in WEEK 3
		return graph.size();
	}

	/**
	 * Return the intersections, which are the vertices in this graph.
	 * 
	 * @return The vertices in this graph as GeographicPoints
	 */
	public Set<GeographicPoint> getVertices() {
		// TODO: Implement this method in WEEK 3
		return new HashSet<GeographicPoint>(graph.keySet());
	}

	/**
	 * Get the number of road segments in the graph
	 * 
	 * @return The number of edges in the graph.
	 */
	public int getNumEdges() {
		// TODO: Implement this method in WEEK 3
		return numEdges;
	}

	/**
	 * Add a node corresponding to an intersection at a Geographic Point If the
	 * location is already in the graph or null, this method does not change the
	 * graph.
	 * 
	 * @param location
	 *            The location of the intersection
	 * @return true if a node was added, false if it was not (the node was
	 *         already in the graph, or the parameter is null).
	 */
	public boolean addVertex(GeographicPoint location) {
		// TODO: Implement this method in WEEK 3
		if (location == null || graph.containsKey(location)) {
			return false;
		}
		graph.put(location, new LinkedList<>());
		return true;
	}

	/**
	 * Adds a directed edge to the graph from pt1 to pt2. Precondition: Both
	 * GeographicPoints have already been added to the graph
	 * 
	 * @param from
	 *            The starting point of the edge
	 * @param to
	 *            The ending point of the edge
	 * @param roadName
	 *            The name of the road
	 * @param roadType
	 *            The type of the road
	 * @param length
	 *            The length of the road, in km
	 * @throws IllegalArgumentException
	 *             If the points have not already been added as nodes to the
	 *             graph, if any of the arguments is null, or if the length is
	 *             less than 0.
	 */
	public void addEdge(GeographicPoint from, GeographicPoint to, String roadName, String roadType, double length)
			throws IllegalArgumentException {

		// System.out.println("edge -> from: (" + from + ") to: (" + to + ")
		// length:" + length);

		// TODO: Implement this method in WEEK 3
		if (!graph.containsKey(from) || !graph.containsKey(to)) {
			throw new IllegalArgumentException("Either from or to is not in graph already");
		}

		if (roadName == null || roadType == null || length < 0) {
			throw new IllegalArgumentException("Road name and type can't be null, " + "and length must be positive");
		}

		List<StreetSegment> streetSegmentList = graph.get(from);

		if (streetSegmentList == null) {
			throw new IllegalArgumentException("streetSegmentList is null, shouldn't be null");
		}

		streetSegmentList.add(new StreetSegment(from, to, roadName, roadType, length));

		numEdges++;
	}

	/**
	 * Find the path from start to goal using breadth first search
	 * 
	 * @param start
	 *            The starting location
	 * @param goal
	 *            The goal location
	 * @return The list of intersections that form the shortest (unweighted)
	 *         path from start to goal (including both start and goal).
	 */
	public List<GeographicPoint> bfs(GeographicPoint start, GeographicPoint goal) {
		// Dummy variable for calling the search algorithms
		Consumer<GeographicPoint> temp = (x) -> {
		};
		return bfs(start, goal, temp);
	}

	/**
	 * Find the path from start to goal using breadth first search.
	 * 
	 * The shortest path BPS solution is a bit inefficient due to the following
	 * reasons: 1) It only considers the number of hops and not the distance of
	 * each path (edge) 2) It explores all the paths of each of the vertexes and
	 * therefore it will take a while. 3) It will blindly explore paths that
	 * will not lead to the destination (not too smart)
	 * 
	 * @param start
	 *            The starting location
	 * @param goal
	 *            The goal location
	 * @param nodeSearched
	 *            A hook for visualization. See assignment instructions for how
	 *            to use it.
	 * @return The list of intersections that form the shortest (unweighted)
	 *         path from start to goal (including both start and goal).
	 */
	public List<GeographicPoint> bfs(GeographicPoint start, GeographicPoint goal,
			Consumer<GeographicPoint> nodeSearched) {
		// TODO: Implement this method in WEEK 3

		if (start == null || goal == null) {
			throw new IllegalArgumentException("start and goal can't be null");
		}

		// Hook for visualization. See writeup.
		// nodeSearched.accept(next.getLocation());

		// using a queue
		Queue<GeographicPoint> queue = new LinkedList<>();
		queue.offer(start);

		Map<GeographicPoint, List<GeographicPoint>> pathMap = new HashMap<>();

		boolean foundStartPoint = false;
		boolean foundIt = false;
		List<GeographicPoint> shorttestPathList = null;
		while (!queue.isEmpty()) {
			GeographicPoint point = queue.poll();
			nodeSearched.accept(point);

			System.out.println("Processing point: " + point);
			nodeSearched.accept(point);

			if (point.equals(goal)) {
				foundIt = true;
				shorttestPathList = pathMap.get(point);
				// shorttestPathList.add(goal);
				break;
			}

			if (point.equals(start)) {
				foundStartPoint = true;
				List<GeographicPoint> pathListSoFar = new LinkedList<>();
				pathListSoFar.add(start);
				pathMap.put(start, pathListSoFar);
			}

			// exploring edges from the current point
			// we do this regardless whether we've found the start point for not
			List<StreetSegment> edgeList = graph.get(point);
			if (edgeList != null) {
				List<GeographicPoint> pathListSoFar = pathMap.get(point);
				for (StreetSegment segment : edgeList) {
					GeographicPoint newVertex = segment.getPoint2();
					// check to see if we have seen this vertex before
					if (!pathMap.containsKey(newVertex)) {
						queue.add(newVertex);

						if (foundStartPoint) {
							// clone it
							List<GeographicPoint> newPathListSoFar = new LinkedList(pathListSoFar);
							newPathListSoFar.add(newVertex);
							pathMap.put(newVertex, newPathListSoFar);
						}
					} else {
						System.out.println("Already visited vertext: " + newVertex);
					}
				}
			}
		}
		return shorttestPathList;
	}

	/**
	 * Find the path from start to goal using Dijkstra's algorithm.
	 * 
	 * 
	 * @param start
	 *            The starting location
	 * @param goal
	 *            The goal location
	 * @return The list of intersections that form the shortest path from start
	 *         to goal (including both start and goal).
	 */
	public List<GeographicPoint> dijkstra(GeographicPoint start, GeographicPoint goal) {
		
		// Dummy variable for calling the search algorithms
		// You do not need to change this method.
		Consumer<GeographicPoint> temp = (x) -> {
		};
		return dijkstra(start, goal, temp);
	}

	/**
	 * Find the path from start to goal using Dijkstra's algorithm.
	 * 
	 * Optimizations: 1) Dijkstra's algorithm uses the weight of reach edge to
	 * intelligent determine which path to explore 2) It still uses BFS, but in
	 * a more optimized way by selecting the shortest path to explore
	 * incrementally
	 * 
	 * 
	 * @param start
	 *            The starting location
	 * @param goal
	 *            The goal location
	 * @param nodeSearched
	 *            A hook for visualization. See assignment instructions for how
	 *            to use it.
	 * @return The list of intersections that form the shortest path from start
	 *         to goal (including both start and goal).
	 */
	public List<GeographicPoint> dijkstra(GeographicPoint start, GeographicPoint goal,
			Consumer<GeographicPoint> nodeSearched) {
		System.out.println("******** dijkstra *****");
		// TODO: Implement this method in WEEK 4

		// Hook for visualization. See write up.
		// nodeSearched.accept(next.getLocation());

		// instead of a regular queue we will need priority queue that sorts
		// items by weight in ascending
		// order. Each item in the priority queue contains a GeographicPoint and
		// associated distance
		PriorityQueue<WeightGeographicPoint> pq = new PriorityQueue<>(new Comparator<WeightGeographicPoint>() {
			@Override
			public int compare(WeightGeographicPoint o1, WeightGeographicPoint o2) {
				double diff = o1.getWeight() - o2.getWeight();
				if (diff > 0) {
					return 1;
				} else if ( diff < 0) {
					return -1;
				} else {
					return 0;
				}
				// the line below will cause overflow problem
				//return (int) (o1.getWeight() - o2.getWeight());
			}
		});

		// we also need a map to easily find out which vertexes have already
		// visited.
		// In addition this map contains the associated weight from starting
		// point and the weight
		// will adjusted as additional paths are explored
		Set<GeographicPoint> visitedVertexSet = new HashSet<>();

		// We also need a map track the shortest distance from start to every
		// vertex
		Map<GeographicPoint, Double> vertexDistanceMap = new HashMap<>();

		// we also need a parent map so we can trace back the shortest path
		Map<GeographicPoint, GeographicPoint> parentMap = new HashMap<>();

		// start
		pq.add(WeightGeographicPoint.create(start, 0));

		vertexDistanceMap.put(start, 0.0);

		boolean foundIt = false;
		int visitCount = 0;
		while (!pq.isEmpty()) {
			WeightGeographicPoint currWeightGP = pq.poll();
			visitCount++;
			
			//System.out.println("visiting: " + currWeightGP);
			GeographicPoint curr = currWeightGP.getGeographicPoint();
			double currentGPWeight = currWeightGP.getWeight();

			nodeSearched.accept(curr);

			// check see if we already have already seen this vertex before
			if (!visitedVertexSet.contains(curr)) {
				// very important to mark it as visited
				visitedVertexSet.add(curr);

				if (curr.equals(goal)) {
					foundIt = true;
					break;
				}

				// for each of the edges from curr and not visited
				List<StreetSegment> edgeList = graph.get(curr);
				if (edgeList != null) {
					for (StreetSegment segment : edgeList) {
						GeographicPoint nextVertex = segment.getPoint2();
						// again, only if a new reachable hasn't been visited
						if (!visitedVertexSet.contains(nextVertex)) {
							// updating the weight if shorter
							double nextVertexCurrentWeight = getVertexCurrentWeight(nextVertex, vertexDistanceMap);
							double newWeight = currentGPWeight + segment.getLength();

							// perform the following only if the new weight is
							// less than previous weight
							if (newWeight < nextVertexCurrentWeight) {
								vertexDistanceMap.put(nextVertex, newWeight);

								// update parent map
								parentMap.put(nextVertex, curr);

								// add new vertexes into priority queue
								pq.add(WeightGeographicPoint.create(nextVertex, newWeight));
							}
						}
					}
				}
			}
		}

		List<GeographicPoint> shortestPathList = new LinkedList<>();
		if (foundIt) {
			// build the path list based on the parent map to return
			Stack<GeographicPoint> pathStack = new Stack<>();
			// walking backward from the goal vertex
			GeographicPoint startVertexToBacktrack = goal;
			while (startVertexToBacktrack != null) {
				pathStack.push(startVertexToBacktrack);
				startVertexToBacktrack = parentMap.get(startVertexToBacktrack);
			}

			// now convert the stack to list with correct order of starting from
			// start node
			while (!pathStack.isEmpty()) {
				shortestPathList.add(pathStack.pop());
			}
		}

		System.out.println("dijkstra visitCount: " + visitCount);
		return shortestPathList;
	}

	private double getVertexCurrentWeight(GeographicPoint nextVertex, Map<GeographicPoint, Double> vertexDistanceMap) {
		double nextVertexCurrentWeight = Double.MAX_VALUE;

		if (vertexDistanceMap.containsKey(nextVertex)) {
			nextVertexCurrentWeight = vertexDistanceMap.get(nextVertex);
		}

		return nextVertexCurrentWeight;
	}

	/**
	 * Find the path from start to goal using A-Star search
	 * 
	 * @param start
	 *            The starting location
	 * @param goal
	 *            The goal location
	 * @return The list of intersections that form the shortest path from start
	 *         to goal (including both start and goal).
	 */
	public List<GeographicPoint> aStarSearch(GeographicPoint start, GeographicPoint goal) {
		// Dummy variable for calling the search algorithms
		Consumer<GeographicPoint> temp = (x) -> {
		};
		return aStarSearch(start, goal, temp);
	}

	/**
	 * Find the path from start to goal using A-Star search. This is an
	 * improvement on the Dijkstra's algorithm by add an additional heuristics
	 * in the prioritizing which node should be visited next.
	 * 
	 * Observe that Dijkstra's doesn't have an eye toward the goal, it just
	 * blindly explores nodes that are not helping in moving to the direction of
	 * goal. This is due to single heuristic that indicates which node is the
	 * closest to explore next, even though it may be going the opposite
	 * direction of the goal.
	 * 
	 * A star algorithm basically adds an additional heuristic help directing
	 * the path toward the end goal and this will reduce the number of nodes the
	 * algorithm needs to visit.
	 * 
	 * Don't forget that A-star algorithm is built upon the Dijkstra's
	 * algorithm. It is not a brand new one.
	 * 
	 * One difference between Dijkstra's algorithm and A-start algorithm is: *
	 * Dijkstra’s Algorithm can find paths to all locations * A* finds paths to
	 * one location
	 * 
	 * With Breadth First Search and Dijkstra’s Algorithm, the frontier expands
	 * in all directions. A* is the best of both worlds. As long as the
	 * heuristic does not overestimate distances
	 * 
	 * Resources:
	 * http://theory.stanford.edu/~amitp/GameProgramming/AStarComparison.html
	 * http://www.redblobgames.com/pathfinding/a-star/introduction.html
	 * 
	 * @param start
	 *            The starting location
	 * @param goal
	 *            The goal location
	 * @param nodeSearched
	 *            A hook for visualization. See assignment instructions for how
	 *            to use it.
	 * @return The list of intersections that form the shortest path from start
	 *         to goal (including both start and goal).
	 */
	public List<GeographicPoint> aStarSearch(GeographicPoint start, GeographicPoint goal,
			Consumer<GeographicPoint> nodeSearched) {
		System.out.println("******** aStarSearch *****");
		// TODO: Implement this method in WEEK 4

		// Hook for visualization. See writeup.
		// nodeSearched.accept(next.getLocation());

		// instead of a regular queue we will need priority queue that sorts
		// items by weight in ascending
		// order. Each item in the priority queue contains a GeographicPoint and
		// associated distance
		PriorityQueue<WeightGeographicPoint> pq = new PriorityQueue<>(new Comparator<WeightGeographicPoint>() {
			@Override
			public int compare(WeightGeographicPoint o1, WeightGeographicPoint o2) {
				double combinedHeuristic1 = o1.getCombinedWeight();
				double combinedHeuristic2 = o2.getCombinedWeight();
				
				double difference = combinedHeuristic1 - combinedHeuristic2;
				
				if (difference < 0.0) {
					return -1;
				} else if (difference > 0.0) {
					return 1;
				} else {
					return 0;
				}
			}
		});

		// we also need a map to easily find out which vertexes have already
		// visited.
		// In addition this map contains the associated weight from starting
		// point and the weight
		// will adjusted as additional paths are explored
		Set<GeographicPoint> visitedVertexSet = new HashSet<>();

		// We also need a map track the shortest distance from start to every
		// vertex
		Map<GeographicPoint, Double> vertexDistanceMap = new HashMap<>();

		// we also need a parent map so we can trace back the shortest path
		Map<GeographicPoint, GeographicPoint> parentMap = new HashMap<>();

		// start
		pq.add(WeightGeographicPoint.create(start, 0, 0));

		vertexDistanceMap.put(start, 0.0);

		boolean foundIt = false;
		int visitCount = 0;
		while (!pq.isEmpty()) {
			WeightGeographicPoint currWeightGP = pq.poll();
			visitCount++;
			//System.out.println("visiting: " + currWeightGP);

			GeographicPoint curr = currWeightGP.getGeographicPoint();
			double currentGPWeight = currWeightGP.getWeight();

			nodeSearched.accept(curr);

			// check see if we already have already seen this vertex before
			if (!visitedVertexSet.contains(curr)) {
				// very important to mark it as visited
				visitedVertexSet.add(curr);

				if (curr.equals(goal)) {
					foundIt = true;
					break;
				}

				// for each of the edges from curr and not visited
				List<StreetSegment> edgeList = graph.get(curr);
				if (edgeList != null) {
					for (StreetSegment segment : edgeList) {
						GeographicPoint nextVertex = segment.getPoint2();
						// again, only if a new reachable hasn't been visited
						if (!visitedVertexSet.contains(nextVertex)) {
							// updating the weight if shorter
							double nextVertexCurrentWeight = getVertexCurrentWeight(nextVertex, vertexDistanceMap);
							double newWeight = currentGPWeight + segment.getLength();

							// perform the following only if the new weight is
							// less than previous weight
							if (newWeight < nextVertexCurrentWeight) {
								vertexDistanceMap.put(nextVertex, newWeight);

								// update parent map
								parentMap.put(nextVertex, curr);

								// add new vertexes into priority queue
								double predictedDistance = nextVertex.distance(goal);
								pq.add(WeightGeographicPoint.create(nextVertex, newWeight, predictedDistance));
							}
						}
					}
				}
			}
		}

		List<GeographicPoint> shortestPathList = new LinkedList<>();
		if (foundIt) {
			// build the path list based on the parent map to return
			Stack<GeographicPoint> pathStack = new Stack<>();
			// walking backward from the goal vertex
			GeographicPoint startVertexToBacktrack = goal;
			while (startVertexToBacktrack != null) {
				pathStack.push(startVertexToBacktrack);
				startVertexToBacktrack = parentMap.get(startVertexToBacktrack);
			}

			// now convert the stack to list with correct order of starting from
			// start node
			while (!pathStack.isEmpty()) {
				shortestPathList.add(pathStack.pop());
			}
		}

		System.out.println("astar visitCount: " + visitCount);
		return shortestPathList;
	}

	public static void main(String[] args) {
		System.out.print("Making a new map...");
		MapGraph firstMap = new MapGraph();
		System.out.print("DONE. \nLoading the map...");
		GraphLoader.loadRoadMap("data/testdata/simpletest.map", firstMap);
		System.out.println("DONE.");

		// You can use this method for testing.

		/*
		 * Here are some test cases you should try before you attempt the Week 3
		 * End of Week Quiz, EVEN IF you score 100% on the programming
		 * assignment.
		 */

		MapGraph simpleTestMap = new MapGraph();
		GraphLoader.loadRoadMap("data/testdata/simpletest.map", simpleTestMap);

		GeographicPoint testStart = new GeographicPoint(1.0, 1.0);
		GeographicPoint testEnd = new GeographicPoint(8.0, -1.0);

		System.out.println("Test 1 using simpletest: Dijkstra should be 9 and AStar should be 5");
		List<GeographicPoint> testroute = simpleTestMap.bfs(testStart, testEnd);
		System.out.println("bfs: testroute: " + testroute);
		
		System.out.println("Dijkstra: testroute: " + simpleTestMap.dijkstra(testStart, testEnd));
		System.out.println("aStarSearch: testroute: " + simpleTestMap.aStarSearch(testStart, testEnd));

		// List<GeographicPoint> testroute =
		// simpleTestMap.dijkstra(testStart,testEnd);
		// List<GeographicPoint> testroute2 =
		// simpleTestMap.aStarSearch(testStart,testEnd);

		//test1();
		
		quiz();
		
		/*
		 * 
		 * MapGraph testMap = new MapGraph();
		 * GraphLoader.loadRoadMap("data/maps/utc.map", testMap);
		 * 
		 * // A very simple test using real data testStart = new
		 * GeographicPoint(32.869423, -117.220917); testEnd = new
		 * GeographicPoint(32.869255, -117.216927); System.out.
		 * println("Test 2 using utc: Dijkstra should be 13 and AStar should be 5"
		 * ); testroute = testMap.dijkstra(testStart,testEnd); testroute2 =
		 * testMap.aStarSearch(testStart,testEnd);
		 * 
		 * 
		 * // A slightly more complex test using real data testStart = new
		 * GeographicPoint(32.8674388, -117.2190213); testEnd = new
		 * GeographicPoint(32.8697828, -117.2244506); System.out.
		 * println("Test 3 using utc: Dijkstra should be 37 and AStar should be 10"
		 * ); testroute = testMap.dijkstra(testStart,testEnd); testroute2 =
		 * testMap.aStarSearch(testStart,testEnd);
		 */

		/* Use this code in Week 3 End of Week Quiz */
		/*
		 * MapGraph theMap = new MapGraph();
		 * System.out.print("DONE. \nLoading the map...");
		 * GraphLoader.loadRoadMap("data/maps/utc.map", theMap);
		 * System.out.println("DONE.");
		 * 
		 * GeographicPoint start = new GeographicPoint(32.8648772,
		 * -117.2254046); GeographicPoint end = new GeographicPoint(32.8660691,
		 * -117.217393);
		 * 
		 * 
		 * List<GeographicPoint> route = theMap.dijkstra(start,end);
		 * List<GeographicPoint> route2 = theMap.aStarSearch(start,end);
		 * 
		 */

	}
	
	private static void test1() {
		 MapGraph testMap = new MapGraph();
		 GraphLoader.loadRoadMap("data/maps/utc.map", testMap);
		  
		  // A very simple test using real data 
		 GeographicPoint testStart = new GeographicPoint(32.869423, -117.220917); 
		 GeographicPoint testEnd = new GeographicPoint(32.869255, -117.216927); 
		 System.out.println("Test 2 using utc: Dijkstra should be 13 and AStar should be 5"); 
		 List<GeographicPoint> testroute = testMap.dijkstra(testStart,testEnd); 
		 List<GeographicPoint> testroute2 = testMap.aStarSearch(testStart,testEnd);
		 
		 System.out.println("======= A very simple test =====");
		 System.out.println("Test 2 using utc: Dijkstra should be 13 and AStar should be 5");
		 System.out.println("There are " + testroute.size() + " nodes for Dijkstra");
		 System.out.println("testroute from dijkstra: " + testroute);
		 System.out.println("There are " + testroute2.size() + " nodes for aStarSearch");
		 System.out.println("testroute from aStarSearch: " + testroute2);
		 
		 System.out.println("");
		 System.out.println("======= A more complex test =====");
		 System.out.println("Test 3 using utc: Dijkstra should be 37 and AStar should be 10");
		 
		 GeographicPoint testStart2 = new GeographicPoint(32.8674388, -117.2190213); 
		 GeographicPoint testEnd2 = new GeographicPoint(32.8697828, -117.2244506); 
		  
		 testroute = testMap.dijkstra(testStart2,testEnd2); 
		 testroute2 = testMap.aStarSearch(testStart2,testEnd2);
		 
		 System.out.println("There are " + testroute.size() + " nodes for Dijkstra");
		 System.out.println("testroute from dijkstra: " + testroute);
		 System.out.println("There are " + testroute2.size() + " nodes for aStarSearch");
		 System.out.println("testroute from aStarSearch: " + testroute2);
		 
	}
	
	private static void quiz() {
		System.out.println("********* auiz **********");
		MapGraph theMap = new MapGraph();
		System.out.print("DONE. \nLoading the map...");
		GraphLoader.loadRoadMap("data/maps/utc.map", theMap);
		System.out.println("DONE.");

		GeographicPoint start = new GeographicPoint(32.8648772, -117.2254046);
		GeographicPoint end = new GeographicPoint(32.8660691, -117.217393);

		List<GeographicPoint> route = theMap.dijkstra(start,end);
		List<GeographicPoint> route2 = theMap.aStarSearch(start,end);
		

	}

}
