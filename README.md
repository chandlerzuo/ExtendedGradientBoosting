# ExtendedGradientBoosting
Extended Gradient Boosting algorithm in c++.

# Description
This package extended the existing implementation of Gradient Boosting algorithm in several ways:
* Use m-of-n rules in buidling individual decision trees;
* Treatment for categorical variables. Instead of creating dummy variables, this implementation seeks the optimal split of all categories in two subsets.
* Mixture modeling. The algorithm learns a mixture of K prediction functions as the final prediction. For individual observations, the algorithm produces the weights across the K clusters, as well as the prediction conditional within each cluster. The final prediction is the weighted sum of individual cluster prediction. 
