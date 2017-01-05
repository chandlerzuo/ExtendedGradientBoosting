# ExtendedGradientBoosting
Extended Gradient Boosting algorithm in c++.
This is a working repository for testing some ideas to improve the existing Gradient Boosting algorithm.
This implementation focus on algorithm improvement and not the computational speed; parallel computation is not supported.

# Description
This package extended the existing implementation of Gradient Boosting algorithm in several ways:

1. **Use m-of-n rules in buidling individual decisition tress.** The decision at each tree branch is based on whether m of n rules hold. Each rule is a binary decision as suggested by [Murphy & Pazzani 91](https://www.cs.rutgers.edu/~pazzani/Publications/mlw91-id2of3.pdf)

2. **Treatment for categorical variables.** Instead of creating dummy variables, this implementation seeks the optimal split of all categories in two subsets. Such splits are learnt to each tree branch. Thus, it creates dynamic binary coding for each category variable. This greatly reduce the feature redundancy, makes feature selection easier, and automatically merge non-informative categories.

3. **Mixture modeling.** The final prediction function is a weighted average of K prediction functions.
  * To do so, the algorithm softly clusters the observations into K groups, and learns a prediction function within each cluster.
  * For prediction, the algorithm calculates the probability of each observation falling into each cluster, and use the probability as the weights to aggregate within-cluster predictions.
  * The cluster probaiblity are soft-max functions of K-1 gradient boosting functions, while the K within-cluster prediction function is also each a gradient boosting function.
  * The algorithm use iteratively gradient descent to simultaneously learn these 2K-1 functions.
