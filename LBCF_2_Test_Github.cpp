#include <Rcpp.h>
using namespace Rcpp;

// Node class definition
class Node {
  
public:
  
  //Attributes
  double mu;
  int variable;
  double split_val;
  LogicalVector observations;
  LogicalVector observations_test;
  bool is_terminal;
  bool in_use;
  
  // Default Constructor
  Node() {
    mu = 0;
    variable = -1;
    split_val = -1;
    is_terminal = false;
    in_use = false;
  }
  
  // Copy constructor for the Node class
  Node(const Node& other) {
    mu = other.mu;
    variable = other.variable;
    split_val = other.split_val;
    observations = clone(other.observations);
    observations_test = clone(other.observations_test);
    is_terminal = other.is_terminal;
    in_use = other.in_use;
  }
  
  
  // Method for updating mu
  void update_mu(double tau, 
                 double tau_mu, 
                 NumericVector y_resid) {
    
    double nj = sum(observations);
    NumericVector node_resid = y_resid[observations];
    double Sj = sum(node_resid);
    mu = R::rnorm((tau * Sj) / (nj * tau + tau_mu), sqrt(1 / (nj * tau + tau_mu)));
  }
  
  
  // Method for updating mu
  void update_tau(double tau, 
                 double tau_tau, 
                 NumericVector y_resid,
                 NumericVector Z) {
    
    NumericVector zy = Z*y_resid;
    NumericVector node_zy_resid = zy[observations];
    NumericVector node_z = Z[observations];
    double nz = sum(node_z);
    double Sj = sum(node_zy_resid);
    mu = R::rnorm((tau * Sj) / (nz * tau + tau_tau), sqrt(1 / (nz * tau + tau_tau)));
  }
  
  
};





// Tree class definition
class Tree {
  
public:
  
  std::vector<Node> node_vector;
  
  // Constructor
  Tree(int num_nodes = 1, int num_obs = 1, int num_obs_test = 1) {
    node_vector.resize(num_nodes);
    node_vector[0].observations = LogicalVector(num_obs, true);
    node_vector[0].observations_test = LogicalVector(num_obs_test, true);
    node_vector[0].in_use=true;
    node_vector[0].is_terminal=true;
  }
  
  // Copy constructor for the Tree class
  Tree(const Tree& other) {
    // Copy each node from the other tree
    for (const Node& node : other.node_vector) {
      node_vector.push_back(Node(node));  // Invoke the Node copy constructor
    }
  }

  
  // Method for updating all terminal nodes
  void update_nodes(double tau, 
                    double tau_mu, 
                    NumericVector y_resid) {
    
    int num_nodes = node_vector.size();
    
    for(int i=0; i<num_nodes; i++)
    {
      if(node_vector[i].is_terminal & node_vector[i].in_use)
      {
        node_vector[i].update_mu(tau, tau_mu, y_resid);
      }
    }
  }
  
  
  // Method for updating all terminal nodes
  void update_nodes_tau(double tau, 
                    double tau_tau, 
                    NumericVector y_resid,
                    NumericVector Z) {
    
    int num_nodes = node_vector.size();
    
    for(int i=0; i<num_nodes; i++)
    {
      if(node_vector[i].is_terminal & node_vector[i].in_use)
      {
        node_vector[i].update_tau(tau, tau_tau, y_resid, Z);
      }
    }
  }

  
  // Method for selecting a terminal node
  int get_terminal_node() {
    
    std::vector<int> valid_indices;
    
    for (int i = 0; i < node_vector.size(); i++) {
      if (node_vector[i].is_terminal & node_vector[i].in_use) {
        valid_indices.push_back(i);
      }
    }
    
    return valid_indices[floor(R::runif(0, valid_indices.size()))];
    
  }
  
  
  // Method for selecting a non terminal node
  int get_non_terminal_node() {
    
    std::vector<int> valid_indices;
    
    for (int i = 0; i < node_vector.size(); i++) {
      if (!node_vector[i].is_terminal & node_vector[i].in_use) {
        valid_indices.push_back(i);
      }
    }
    
    if(valid_indices.size()>0)
    {
      return valid_indices[floor(R::runif(0, valid_indices.size()))];
    }
    else
    {
      return -1;
    }
    
  }
  
  
  // Method for selecting a non terminal node with a parent
  int get_parent_child() {
    
    std::vector<int> valid_indices;
    
    for (int i = 0; i < node_vector.size(); i++) {
      if (!node_vector[i].is_terminal & node_vector[i].in_use & i!=0) {
        valid_indices.push_back(i);
      }
    }
    
    if(valid_indices.size()>0)
    {
      return valid_indices[floor(R::runif(0, valid_indices.size()))];
    }
    else
    {
      return -1;
    }
    
  }
  
  
  // Method for selecting a parent of two terminal nodes
  int get_terminal_parent() {
    
    std::vector<int> valid_indices;
    
    for (int i = 0; i < node_vector.size(); i++) {
      if(node_vector[i].in_use & !node_vector[i].is_terminal)
      {
        if(node_vector[2*i+1].in_use & node_vector[2*i+1].is_terminal & node_vector[2*i+2].in_use & node_vector[2*i+2].is_terminal)
        {
          valid_indices.push_back(i);
        }
      }
    }
    
    if(valid_indices.size()>0)
    {
      return valid_indices[floor(R::runif(0, valid_indices.size()))];
    }
    else
    {
      return -1;
    }
  }
  
  
  
  // Method for growing tree
  void grow(NumericMatrix X, NumericMatrix X_test, int p, int min_nodesize) {
    
    int grow_index = get_terminal_node();
    
    int variable = floor(R::runif(0, p));
    
    node_vector[grow_index].variable = variable;
    
    NumericVector X_col = X(_, variable);
    NumericVector X_col_test = X_test(_, variable);
    
    NumericVector X_col_subset = X_col[node_vector[grow_index].observations];
    NumericVector X_col_subset_test = X_col_test[node_vector[grow_index].observations_test];
    
    NumericVector X_unique = unique(X_col_subset);
    
    double split_val;
    
    if(X_unique.size()>0)
    {
      split_val = sample(X_unique, 1)[0];
    }
    else
    {
      split_val = -1;
    }
    
    node_vector[grow_index].split_val = split_val;
    
    LogicalVector is_less = X_col<=split_val;
    LogicalVector is_less_test = X_col_test<=split_val;
    
    LogicalVector less_subset = node_vector[grow_index].observations & is_less;
    
    LogicalVector more_subset = node_vector[grow_index].observations & !is_less;
    
    int sum_less = sum(less_subset);
    
    int sum_more = sum(more_subset);
    
    if(sum_more>=min_nodesize & sum_less>=min_nodesize)
    {
      if(node_vector.size()<2*grow_index+2+1)
      {
        node_vector.resize(2*grow_index+2+1);
      }
      
      int child_left = 2*grow_index+1;
      
      int child_right = 2*grow_index+2;
      
      node_vector[child_left].observations = node_vector[grow_index].observations & is_less;
      node_vector[child_left].observations_test = node_vector[grow_index].observations_test & is_less_test;
      node_vector[child_left].is_terminal = true;
      node_vector[child_left].in_use = true;
      
      node_vector[child_right].observations = node_vector[grow_index].observations & !is_less;
      node_vector[child_right].observations_test = node_vector[grow_index].observations_test & !is_less_test;
      node_vector[child_right].is_terminal = true;
      node_vector[child_right].in_use = true;
      
      node_vector[grow_index].is_terminal = false;
      node_vector[grow_index].in_use = true;
    }
  }
  
  
  // Method for pruning
  void prune() {
    
    int prune_index = get_terminal_parent();
    
    if(prune_index!=-1)
    {
      node_vector[prune_index*2+1].in_use = false;
      node_vector[prune_index*2+1].is_terminal = false;
      
      node_vector[prune_index*2+2].in_use = false;
      node_vector[prune_index*2+2].is_terminal = false;
      
      node_vector[prune_index].is_terminal = true;
      node_vector[prune_index].in_use = true;
    }
  }
  
  
  
  // Method for updating observations
  void change_update(NumericMatrix X, NumericMatrix X_test) {
    
    int num_nodes = node_vector.size();
    
    for(int i = 0; i < num_nodes; i++)
    {
      if(!node_vector[i].is_terminal & node_vector[i].in_use)
      {
        int child_left = 2*i+1;
        int child_right = 2*i+2;
        
        int variable = node_vector[i].variable;
        double split_val = node_vector[i].split_val;
        
        LogicalVector is_less = X(_, variable)<=split_val;
        LogicalVector is_more = X(_, variable)>split_val;
        
        LogicalVector is_less_test = X_test(_, variable)<=split_val;
        LogicalVector is_more_test = X_test(_, variable)>split_val;
        
        node_vector[child_left].observations = node_vector[i].observations & is_less;
        node_vector[child_right].observations = node_vector[i].observations & is_more;
        
        node_vector[child_left].observations_test = node_vector[i].observations_test & is_less_test;
        node_vector[child_right].observations_test = node_vector[i].observations_test & is_more_test;
      }
    }
  }
  
  
  // Method for changing
  void change(NumericMatrix X, int p) {
    
    int change_index = get_non_terminal_node();
    
    if(change_index!=-1)
    {
      int variable = floor(R::runif(0, p));
      
      node_vector[change_index].variable = variable;
      
      NumericVector X_col = X(_, variable);
      
      X_col = X_col[node_vector[change_index].observations];
      
      NumericVector X_unique = unique(X_col);
      
      if(X_unique.size()>0)
      {
        node_vector[change_index].split_val = sample(X_unique, 1)[0];
      }
      else
      {
        node_vector[change_index].split_val = -1;
      }
    }
  }

  
  
  // Method for swapping
  void swap() {
    
    int swap_index = get_parent_child();
    
    if(swap_index!=-1)
    {
      int parent_index = (swap_index-1)/2;
      
      int parent_variable = node_vector[parent_index].variable;
      double parent_split_val = node_vector[parent_index].split_val;
      
      int child_variable = node_vector[swap_index].variable;
      double child_split_val = node_vector[swap_index].split_val;
      
      node_vector[parent_index].variable = child_variable;
      node_vector[parent_index].split_val = child_split_val;
      
      node_vector[swap_index].variable = parent_variable;
      node_vector[swap_index].split_val = parent_split_val;
    }
  }
  
  
  // Method for checking if any nodes are empty
  bool has_empty_nodes(int min_nodesize) {
    
    int num_nodes = node_vector.size();
    
    for(int i=0; i<num_nodes; i++)
    {
      if(node_vector[i].in_use & node_vector[i].is_terminal)
      {
        if(sum(node_vector[i].observations)<min_nodesize)
        {
          return true;
        } 
      }
    }
    
    return false;
  }

  double log_lik(double tau_mu, 
                 double tau,
                 double alpha,
                 double beta,
                 NumericVector y_resid){
    
    
    double log_lik = 0.0;
    
    for(int i = 0; i < node_vector.size(); i++)
    {
      if(node_vector[i].in_use & node_vector[i].is_terminal)
      {
        double nj = sum(node_vector[i].observations);
        NumericVector node_resid = y_resid[node_vector[i].observations];
        double sum_Rji2 = sum(node_resid*node_resid);
        double Rj_bar = mean(node_resid);
        
        double eq1 = (nj/2.0)*log(tau) + (1.0/2.0)*log(tau_mu/(tau_mu + nj*tau)) - (tau/2.0)*(sum_Rji2 - (tau*(nj*Rj_bar)*(nj*Rj_bar))/(tau_mu + nj*tau));
        double eq4p1 = log(1.0-alpha*pow(1+floor(log2(i + 1)), (-1*beta)));
        
        log_lik += eq1 + eq4p1;
      }
      else if(node_vector[i].in_use & !node_vector[i].is_terminal)
      {
        double eq4p2 = log(alpha)-beta*log(1+floor(log2(i + 1)));
        
        log_lik += eq4p2;
      }
    }
    
    return log_lik;
  }
  
  
  double log_lik_tau(double tau_tau, 
                 double tau,
                 double alpha_tau,
                 double beta_tau,
                 NumericVector y_resid,
                 NumericVector Z){
    
    
    double log_lik = 0.0;
    
    for(int i = 0; i < node_vector.size(); i++)
    {
      if(node_vector[i].in_use & node_vector[i].is_terminal)
      {
        double nj = sum(node_vector[i].observations);
        NumericVector node_resid = y_resid[node_vector[i].observations];
        NumericVector zy = Z*y_resid;
        NumericVector node_zy_resid= zy[node_vector[i].observations];
        NumericVector node_z = Z[node_vector[i].observations];
        double nz = sum(node_z);
        double sum_Rji2 = sum(node_resid*node_resid);
        
        double eq1 = (nj/2.0)*log(tau) + (1.0/2.0)*log(tau_tau/(tau_tau + nz*tau)) - (tau/2.0)*(sum_Rji2 - (tau*(sum(node_zy_resid))*(sum(node_zy_resid)))/(tau_tau + nz*tau));
        double eq4p1 = log(1.0-alpha_tau*pow(1+floor(log2(i + 1)), (-1*beta_tau)));
        
        log_lik += eq1 + eq4p1;
      }
      else if(node_vector[i].in_use & !node_vector[i].is_terminal)
      {
        double eq4p2 = log(alpha_tau)-beta_tau*log(1+floor(log2(i + 1)));
        
        log_lik += eq4p2;
      }
    }
    
    return log_lik;
  }
  
  


  
  // Method for getting predictions from tree
  NumericVector get_predictions() {
    int num_obs = node_vector[0].observations.size();
    int num_nodes = node_vector.size();
    
    NumericVector predictions(num_obs, 0.0);
    
    for(int i = 0; i < num_nodes; i++) {
      for(int j = 0; j < num_obs; j++) {
        if(node_vector[i].is_terminal & node_vector[i].in_use & node_vector[i].observations[j]) {
          predictions[j] = node_vector[i].mu;
        }
      }
    }
    
    return predictions;
  }
  
  
  // Method for getting predictions from tree
  NumericVector get_predictions_test() {
    int num_obs = node_vector[0].observations_test.size();
    int num_nodes = node_vector.size();
    
    NumericVector predictions_test(num_obs, 0.0);
    
    for(int i = 0; i < num_nodes; i++) {
      for(int j = 0; j < num_obs; j++) {
        if(node_vector[i].is_terminal & node_vector[i].in_use & node_vector[i].observations_test[j]) {
          predictions_test[j] = node_vector[i].mu;
        }
      }
    }
    
    return predictions_test;
  }
  
};


// Forest class definition
class Forest {
  
public:
  
  std::vector<Tree> tree_vector;
  
  // Constructor
  Forest(int num_trees=1, int num_nodes = 1, int num_obs=1, int num_obs_test = 1) {
    
    tree_vector.resize(num_trees);
    
    for(int i=0; i<num_trees; i++)
    {
      tree_vector[i] = Tree(num_nodes, num_obs, num_obs_test);
    }
  }
};


NumericVector rowSumsWithoutColumn(NumericMatrix mat, int columnToRemove) {
  int numRows = mat.nrow();
  int numCols = mat.ncol();
  
  // Create a vector to store the row sums
  NumericVector rowSums(numRows, 0.0);
  
  // Iterate over each row
  for (int i = 0; i < numRows; i++) {
    // Iterate over each column, excluding the column to remove
    for (int j = 0; j < numCols; j++) {
      if (j != columnToRemove) {
        // Add the element to the row sum
        rowSums[i] += mat(i, j);
      }
    }
  }
  
  return rowSums;
}
 
 

double sample_tau(double n,
                  double nu,
                  NumericVector y,
                  NumericVector preds,
                  double lambda) {
  
  double shape = (n+nu)/2.0;
  double S = sum(pow((y-preds), 2));
  double rate = (S+nu*lambda)/2.0;
  double scale = 1.0/rate;
  
  return R::rgamma(shape, scale);
}


// [[Rcpp::export]]
List fast_lbcf2(NumericMatrix X_mu,
               NumericVector y,
               NumericVector Z_growth,
               NumericVector Z_treat,
               NumericMatrix X_growth,
               NumericMatrix X_treat,
               NumericMatrix X_mu_test,
               NumericMatrix X_growth_test,
               NumericMatrix X_treat_test,
               double alpha_mu,
               double beta_mu,
               double alpha_growth,
               double beta_growth,
               double alpha_treat,
               double beta_treat,
               double tau_mu,
               double tau_growth,
               double tau_treat,
               double nu,
               double lambda,
               int n_iter,
               int n_tree_mu,
               int n_tree_growth,
               int n_tree_treat,
               int min_nodesize)
{
  //set tau
  double tau = 1;
  
  //normalise y
  double y_mean = mean(y);
  double y_sd = sd(y);
  NumericVector y_scaled = (y-y_mean)/y_sd;
  
  //get number of variables p, and rows n
  int n = y_scaled.size();
  int n_test = X_mu_test.nrow();
  int p_mu = X_mu.ncol();
  int p_growth = X_growth.ncol();
  int p_treat = X_treat.ncol();
  
  //For holding tree predictions at each iteration
  NumericMatrix tree_preds_mu(n, n_tree_mu);
  NumericMatrix tree_preds_growth(n, n_tree_growth);
  NumericMatrix tree_preds_treat(n, n_tree_treat);
  
  NumericMatrix tree_preds_mu_test(n_test, n_tree_mu);
  NumericMatrix tree_preds_growth_test(n_test, n_tree_growth);
  NumericMatrix tree_preds_treat_test(n_test, n_tree_treat);
  
  //For holding overall predictions from each iteration
  NumericMatrix preds_mat_mu(n, n_iter);
  NumericMatrix preds_mat_growth(n, n_iter);
  NumericMatrix preds_mat_treat(n, n_iter);
  
  NumericMatrix preds_mat_mu_test(n_test, n_iter);
  NumericMatrix preds_mat_growth_test(n_test, n_iter);
  NumericMatrix preds_mat_treat_test(n_test, n_iter);
  
  //For holding tau from each iteration
  NumericVector taus(n_iter);
  
  StringVector choices = {"Grow", "Prune", "Change", "Swap"};
  
  Forest mu_forest(n_tree_mu, 1, n, n_test);
  
  Forest growth_forest(n_tree_growth, 1, n, n_test);
  
  Forest treat_forest(n_tree_treat, 1, n, n_test);
  
  for(int iter = 0; iter < n_iter; iter++)
  {
    for(int tree_num = 0; tree_num < n_tree_mu; tree_num++)
    {
      NumericVector y_resid = y_scaled-rowSumsWithoutColumn(tree_preds_mu, tree_num)-Z_growth*rowSumsWithoutColumn(tree_preds_growth, -1)-Z_treat*rowSumsWithoutColumn(tree_preds_treat, -1);
      
      String choice = sample(choices, 1)[0];
      
      Tree proposal_tree = Tree(mu_forest.tree_vector[tree_num]);
      
      if(choice == "Grow")
      {
        proposal_tree.grow(X_mu, X_mu_test, p_mu, min_nodesize);
      }
      
      
      if(choice == "Prune")
      {
        proposal_tree.prune();
      }
      
      
      if(choice == "Change")
      {
        proposal_tree.change(X_mu, p_mu);
        proposal_tree.change_update(X_mu, X_mu_test);
      }
      
      
      if(choice == "Swap")
      {
        proposal_tree.swap();
        proposal_tree.change_update(X_mu, X_mu_test);
      }
      
      
      if(!proposal_tree.has_empty_nodes(min_nodesize))
      {
        double lnew = proposal_tree.log_lik(tau_mu, 
                                            tau,
                                            alpha_mu,
                                            beta_mu,
                                            y_resid);
          
        double lold = mu_forest.tree_vector[tree_num].log_lik(tau_mu, 
                                                                tau,
                                                                alpha_mu,
                                                                beta_mu,
                                                                y_resid);
          
        double a = exp(lnew-lold);
        if(a > R::runif(0, 1))
        {
          mu_forest.tree_vector[tree_num] = Tree(proposal_tree);
        }
      }
      
      mu_forest.tree_vector[tree_num].update_nodes(tau, tau_mu, y_resid);
      
      NumericVector tree_preds_from_iter_mu = mu_forest.tree_vector[tree_num].get_predictions();
      NumericVector tree_preds_from_iter_mu_test = mu_forest.tree_vector[tree_num].get_predictions_test();
      
      for(int i=0; i<n; i++)
      {
        tree_preds_mu(i, tree_num) = tree_preds_from_iter_mu[i];
      }
      for(int i=0; i<n_test; i++)
      {
        tree_preds_mu_test(i, tree_num) = tree_preds_from_iter_mu_test[i];
      }
    
    }
    
    for(int tree_num = 0; tree_num < n_tree_growth; tree_num++)
    {
      NumericVector y_resid = y_scaled-rowSumsWithoutColumn(tree_preds_mu, -1)-Z_growth*rowSumsWithoutColumn(tree_preds_growth, tree_num)-Z_treat*rowSumsWithoutColumn(tree_preds_treat, -1);
      
      String choice = sample(choices, 1)[0];
      
      Tree proposal_tree = Tree(growth_forest.tree_vector[tree_num]);
      
      if(choice == "Grow")
      {
        proposal_tree.grow(X_growth, X_growth_test, p_growth, min_nodesize);
      }
      
      
      if(choice == "Prune")
      {
        proposal_tree.prune();
      }
      
      
      if(choice == "Change")
      {
        proposal_tree.change(X_growth, p_growth);
        proposal_tree.change_update(X_growth, X_growth_test);
      }
      
      
      if(choice == "Swap")
      {
        proposal_tree.swap();
        proposal_tree.change_update(X_growth, X_growth_test);
      }
      
      
      if(!proposal_tree.has_empty_nodes(min_nodesize))
      {
        double lnew = proposal_tree.log_lik_tau(tau_growth, 
                                            tau,
                                            alpha_growth,
                                            beta_growth,
                                            y_resid,
                                            Z_growth);
        
        double lold = growth_forest.tree_vector[tree_num].log_lik_tau(tau_growth, 
                                                                tau,
                                                                alpha_growth,
                                                                beta_growth,
                                                                y_resid,
                                                                Z_growth);
        
        double a = exp(lnew-lold);
        if(a > R::runif(0, 1))
        {
          growth_forest.tree_vector[tree_num] = Tree(proposal_tree);
        }
      }
      
      growth_forest.tree_vector[tree_num].update_nodes_tau(tau, tau_growth, y_resid, Z_growth);
      
      NumericVector tree_preds_from_iter_growth = growth_forest.tree_vector[tree_num].get_predictions();
      NumericVector tree_preds_from_iter_growth_test = growth_forest.tree_vector[tree_num].get_predictions_test();
      
      for(int i=0; i<n; i++)
      {
        tree_preds_growth(i, tree_num) = tree_preds_from_iter_growth[i];
      }
      for(int i=0; i<n_test; i++)
      {
        tree_preds_growth_test(i, tree_num) = tree_preds_from_iter_growth_test[i];
      }
      
    }
    
    
    for(int tree_num = 0; tree_num < n_tree_treat; tree_num++)
    {
      NumericVector y_resid = y_scaled-rowSumsWithoutColumn(tree_preds_mu, -1)-Z_growth*rowSumsWithoutColumn(tree_preds_growth, -1)-Z_treat*rowSumsWithoutColumn(tree_preds_treat, tree_num);
      
      String choice = sample(choices, 1)[0];
      
      Tree proposal_tree = Tree(treat_forest.tree_vector[tree_num]);
      
      if(choice == "Grow")
      {
        proposal_tree.grow(X_treat, X_treat_test, p_treat, min_nodesize);
      }
      
      
      if(choice == "Prune")
      {
        proposal_tree.prune();
      }
      
      
      if(choice == "Change")
      {
        proposal_tree.change(X_treat, p_treat);
        proposal_tree.change_update(X_treat, X_treat_test);
      }
      
      
      if(choice == "Swap")
      {
        proposal_tree.swap();
        proposal_tree.change_update(X_treat, X_treat_test);
      }
      
      
      if(!proposal_tree.has_empty_nodes(min_nodesize))
      {
        double lnew = proposal_tree.log_lik_tau(tau_treat, 
                                                tau,
                                                alpha_treat,
                                                beta_treat,
                                                y_resid,
                                                Z_treat);
        
        double lold = treat_forest.tree_vector[tree_num].log_lik_tau(tau_treat, 
                                                                      tau,
                                                                      alpha_treat,
                                                                      beta_treat,
                                                                      y_resid,
                                                                      Z_treat);
        
        double a = exp(lnew-lold);
        if(a > R::runif(0, 1))
        {
          treat_forest.tree_vector[tree_num] = Tree(proposal_tree);
        }
      }
      
      treat_forest.tree_vector[tree_num].update_nodes_tau(tau, tau_treat, y_resid, Z_treat);
      
      NumericVector tree_preds_from_iter_treat = treat_forest.tree_vector[tree_num].get_predictions();
      NumericVector tree_preds_from_iter_treat_test = treat_forest.tree_vector[tree_num].get_predictions_test();
      
      for(int i=0; i<n; i++)
      {
        tree_preds_treat(i, tree_num) = tree_preds_from_iter_treat[i];
      }
      for(int i=0; i<n_test; i++)
      {
        tree_preds_treat_test(i, tree_num) = tree_preds_from_iter_treat_test[i];
      }
      
    }
      
    Rcpp::Rcout << "Total of " << iter+1 << " of " << n_iter << " iterations completed! " << "(" << (float)(iter+1)/(float)n_iter*100 << "%)             " << "\r";
    Rcpp::Rcout.flush();
    
    NumericVector iter_preds_mu = rowSumsWithoutColumn(tree_preds_mu, -1);
    NumericVector iter_preds_growth = rowSumsWithoutColumn(tree_preds_growth, -1);
    NumericVector iter_preds_treat = rowSumsWithoutColumn(tree_preds_treat, -1);
    
    NumericVector iter_preds_mu_test = rowSumsWithoutColumn(tree_preds_mu_test, -1);
    NumericVector iter_preds_growth_test = rowSumsWithoutColumn(tree_preds_growth_test, -1);
    NumericVector iter_preds_treat_test = rowSumsWithoutColumn(tree_preds_treat_test, -1);
    
    for(int i=0; i<n; i++)
    {
      preds_mat_mu(i, iter) = iter_preds_mu[i];
      preds_mat_growth(i, iter) = iter_preds_growth[i];
      preds_mat_treat(i, iter) = iter_preds_treat[i];
    }
    for(int i=0; i<n_test; i++)
    {
      preds_mat_mu_test(i, iter) = iter_preds_mu_test[i];
      preds_mat_growth_test(i, iter) = iter_preds_growth_test[i];
      preds_mat_treat_test(i, iter) = iter_preds_treat_test[i];
    }
    
    tau=sample_tau(n, nu, y_scaled, iter_preds_mu+Z_growth*iter_preds_growth+Z_treat*iter_preds_treat, lambda);
    
    taus[iter] = tau;
  }
  
  Rcpp::Rcout << "";
  
  
  return List::create(
    Named("predictions_mu") = (preds_mat_mu*y_sd)+y_mean,
    Named("predictions_growth") = preds_mat_growth*y_sd,
    Named("predictions_treat") = preds_mat_treat*y_sd,
    Named("predictions_mu_test") = (preds_mat_mu_test*y_sd)+y_mean,
    Named("predictions_growth_test") = preds_mat_growth_test*y_sd,
    Named("predictions_treat_test") = preds_mat_treat_test*y_sd,
    Named("taus") = taus/(pow(y_sd, 2)),
    Named("sigmas") = y_sd/sqrt(taus)
  );
}
 

