# Master-Assisted Distributed Uplink Operation for Cell-Free Massive MIMO Networks (Accepted Paper ICASSP 2026)
**Authros**: Andreas Angelou, Pourya Behmandpoor, Marc Moonen

The paper link: [https://ieeexplore.ieee.org/abstract/document/11464100](https://ieeexplore.ieee.org/abstract/document/11464100)

**Brief Introduction**: Master-Assisted Distributed Uplink Operation (MADUO) is a novel uplink operation for cell-free massive MIMO networks. Unlike centralized operation, where APs forward their raw antenna signals to the CPU, and distributed operation, where APs send local signal estimates to the CPU, MADUO assigns each UE a set of additional serving APs (ASAPs) and a master AP (MAP). Each ASAP computes a local estimate of the UE’s signal and sends it to the MAP. The MAP then combines the received signal estimates from the ASAPs with its own antenna signals to obtain the final signal estimate.

**Acknowledgments**: Parts of this code are based on the publicly available implementation at [https://github.com/emilbjornson/cell-free-book.git](https://github.com/emilbjornson/cell-free-book.git) and have been further developed by the authors.
