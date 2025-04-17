# A Parallelization of Dijkstra's Shortest Path Algorithm

Implementation of a parallel version of the Dijkstra algorithm by A. Krauser, K. Mehlhorn, U. Meyer, and P. Sanders and subsequent testing. 
n-the number of vertices, m=1.5*n

## What is implemented?

1. Dijkstra's implementation with complexity O(n^2).
2. Implementation of the priority queue.
3. Dijkstra's implementation in O(n*log n) with a written priority queue.
4. Dijkstra's implementation for O(n*log n) with std::priority_queue.
5. Dijkstra's Parallel Algorithm by A. Krauser, K. Mehlhorn, U. Meyer, and P. Sanders.
6. Tests with graphs with sizes from 300 000 to 6 000 000.