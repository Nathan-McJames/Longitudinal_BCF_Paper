//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <RcppDist.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppDist)]]

// Node class definition
class Node {
  
public:
  
  //Attributes
  double mu;
  int variable;
  double split_val;
  LogicalVector observations;
  bool is_terminal;
  bool in_use;
  bool missing_to_left;
  bool use_a_split;
  
  // Default Constructor
  Node() {
    mu = 0;
    variable = -1;
    split_val = -1;
    is_terminal = false;
    in_use = false;
    missing_to_left = true;
    use_a_split = true;
  }
  
  // Copy constructor for the Node class
  Node(const Node& other) {
    mu = other.mu;
    variable = other.variable;
    split_val = other.split_val;
    observations = clone(other.observations);
    is_terminal = other.is_terminal;
    in_use = other.in_use;
    missing_to_left = other.missing_to_left;
    use_a_split = other.use_a_split;
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
  Tree(int num_nodes = 1, int num_obs = 1) {
    node_vector.resize(num_nodes);
    node_vector[0].observations = LogicalVector(num_obs, true);
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
  void grow(NumericMatrix X, int p, int min_nodesize) {
    
    int grow_index = get_terminal_node();
    
    int variable = floor(R::runif(0, p));
    
    node_vector[grow_index].variable = variable;
    
    NumericVector X_col = X(_, variable);
    
    LogicalVector is_missing(node_vector[grow_index].observations.size(), false);
    
    for(int k=0; k<node_vector[grow_index].observations.size(); k++)
    {
      is_missing[k] = !std::isfinite(X_col[k]);
    }
    
    NumericVector X_col_subset = X_col[node_vector[grow_index].observations & !is_missing];
    
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
    
    LogicalVector is_less(node_vector[grow_index].observations.size(), false);
    
    bool missing_to_left = true;
    bool use_a_split = true;
    if(sum(is_missing)==0)
    {
      use_a_split = true;
    }
    else
    {
      if(R::runif(0, 1)>=0.5)
      {
        use_a_split=true;
      }
    }
    
    ////
    //use_a_split = true;
    ////
    
    if(R::runif(0, 1)>0.5)
    {
      missing_to_left = false;
    }
    
    for(int k = 0; k<node_vector[grow_index].observations.size(); k++)
    {
      if(X_col[k]<=split_val)
      {
        is_less[k] = true;
      }
      if(is_missing[k] & missing_to_left)
      {
        is_less[k] = true;
      }
    }
    
    if(!use_a_split)
    {
      for(int k = 0; k<node_vector[grow_index].observations.size(); k++)
      {
        if(is_missing[k])
        {
          is_less[k] = true;
        }
        else{
          is_less[k] = false;
        }
      }
    }
    
    LogicalVector less_subset = node_vector[grow_index].observations & is_less;
    
    LogicalVector more_subset = node_vector[grow_index].observations & !is_less;
    
    int sum_less = sum(less_subset);
    
    int sum_more = sum(more_subset);
    
    node_vector[grow_index].use_a_split = use_a_split;
    node_vector[grow_index].missing_to_left = missing_to_left;
    
    if(sum_more>=min_nodesize & sum_less>=min_nodesize)
    {
      if(node_vector.size()<2*grow_index+2+1)
      {
        node_vector.resize(2*grow_index+2+1);
      }
      
      int child_left = 2*grow_index+1;
      
      int child_right = 2*grow_index+2;
      
      node_vector[child_left].observations = node_vector[grow_index].observations & is_less;
      node_vector[child_left].is_terminal = true;
      node_vector[child_left].in_use = true;
      
      node_vector[child_right].observations = node_vector[grow_index].observations & !is_less;
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
  void change_update(NumericMatrix X) {
    
    int num_nodes = node_vector.size();
    
    for(int i = 0; i < num_nodes; i++)
    {
      if(!node_vector[i].is_terminal & node_vector[i].in_use)
      {
        int child_left = 2*i+1;
        int child_right = 2*i+2;
        
        int variable = node_vector[i].variable;
        double split_val = node_vector[i].split_val;
        
        bool use_a_split = node_vector[i].use_a_split;
        bool missing_to_left = node_vector[i].missing_to_left;
        
        LogicalVector is_less(node_vector[i].observations.size(), false);
        
        LogicalVector is_missing(node_vector[i].observations.size(), false);
        
        for(int k=0; k<node_vector[i].observations.size(); k++)
        {
          is_missing[k] = !std::isfinite(X(k, variable));
        }
        
        for(int k = 0; k<node_vector[i].observations.size(); k++)
        {
          if(X(k, variable)<=split_val)
          {
            is_less[k] = true;
          }
          if(is_missing[k] & missing_to_left)
          {
            is_less[k] = true;
          }
        }
        
        if(!use_a_split)
        {
          for(int k = 0; k<node_vector[i].observations.size(); k++)
          {
            if(is_missing[k])
            {
              is_less[k] = true;
            }
            else{
              is_less[k] = false;
            }
          }
        }
        
        node_vector[child_left].observations = node_vector[i].observations & is_less;
        node_vector[child_right].observations = node_vector[i].observations & !is_less;
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
      
      LogicalVector is_missing(node_vector[change_index].observations.size(), false);
      
      for(int k=0; k<node_vector[change_index].observations.size(); k++)
      {
        is_missing[k] = !std::isfinite(X(k, variable));
      }
      
      NumericVector X_col = X(_, variable);
      
      NumericVector X_col_subset = X_col[node_vector[change_index].observations & !is_missing];
      
      NumericVector X_unique = unique(X_col_subset);
      
      if(X_unique.size()>0)
      {
        node_vector[change_index].split_val = sample(X_unique, 1)[0];
      }
      else
      {
        node_vector[change_index].split_val = -1;
      }
      
      bool missing_to_left = true;
      bool use_a_split = true;
      if(sum(is_missing)==0)
      {
        use_a_split = true;
      }
      else
      {
        if(R::runif(0, 1)>=0.5)
        {
          use_a_split=true;
        }
      }
      
      ////
      //use_a_split = true;
      ////
      
      if(R::runif(0, 1)>0.5)
      {
        missing_to_left = false;
      }
      
      node_vector[change_index].use_a_split = use_a_split;
      node_vector[change_index].missing_to_left = missing_to_left;
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
      bool parent_use_a_split = node_vector[parent_index].use_a_split;
      bool parent_missing_to_left = node_vector[parent_index].missing_to_left;
      
      int child_variable = node_vector[swap_index].variable;
      double child_split_val = node_vector[swap_index].split_val;
      bool child_use_a_split = node_vector[swap_index].use_a_split;
      bool child_missing_to_left = node_vector[swap_index].missing_to_left;
      
      node_vector[parent_index].variable = child_variable;
      node_vector[parent_index].split_val = child_split_val;
      node_vector[parent_index].use_a_split = child_use_a_split;
      node_vector[parent_index].missing_to_left = child_missing_to_left;
      
      node_vector[swap_index].variable = parent_variable;
      node_vector[swap_index].split_val = parent_split_val;
      node_vector[swap_index].use_a_split = parent_use_a_split;
      node_vector[swap_index].missing_to_left = parent_missing_to_left;
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
  
};


// Forest class definition
class Forest {
  
public:
  
  std::vector<Tree> tree_vector;
  
  // Constructor
  Forest(int num_trees=1, int num_nodes = 1, int num_obs=1) {
    
    tree_vector.resize(num_trees);
    
    for(int i=0; i<num_trees; i++)
    {
      tree_vector[i] = Tree(num_nodes, num_obs);
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


//For updating missing z values
double prob_z(double prior_prob,
              double current_precision,
              double y_val,
              double mu_val,
              double growth_val,
              double treat_val)
{
  double b1 = (1-prior_prob)/prior_prob;
  
  double resid_with_z_squared = pow((y_val-mu_val-growth_val-treat_val), 2);
  
  double resid_without_z_squared = pow((y_val-mu_val-growth_val), 2);
  
  double b2 = exp((current_precision/2)*(resid_with_z_squared-resid_without_z_squared));
  
  return 1/(1+b1*b2);
}


// [[Rcpp::export]]
List fast_lbcf2(NumericMatrix X_mu,
               NumericVector y,
               NumericVector Z_growth,
               NumericVector Z_treat,
               NumericMatrix X_growth,
               NumericMatrix X_treat,
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
               int min_nodesize,
               NumericVector p_scores,
               LogicalVector z_missing,
               int n_burn,
               int keep_every)
{
  
  //how many to keep
  int num_keep = 0;
  int num_kept = 0;
  for(int i = 0; i<n_iter; i++)
  {
    if(i>=(n_burn-1) & (i-n_burn) % keep_every==0)
    {
      num_keep++;
    }
  }
  
  
  //set tau
  double tau = 1;
  
  //normalise y
  double y_mean = mean(y);
  double y_sd = sd(y);
  NumericVector y_scaled = (y-y_mean)/y_sd;
  
  //get number of variables p, and rows n
  int n = y_scaled.size();
  int p_mu = X_mu.ncol();
  int p_growth = X_growth.ncol();
  int p_treat = X_treat.ncol();
  
  //For holding tree predictions at each iteration
  NumericMatrix tree_preds_mu(n, n_tree_mu);
  NumericMatrix tree_preds_growth(n, n_tree_growth);
  NumericMatrix tree_preds_treat(n, n_tree_treat);
  
  //For holding overall predictions from each iteration
  NumericMatrix preds_mat_mu(n, num_keep);
  NumericMatrix preds_mat_growth(n, num_keep);
  NumericMatrix preds_mat_treat(n, num_keep);
  
  //For holding tau from each iteration
  NumericVector taus(num_keep);
  
  StringVector choices = {"Grow", "Prune", "Change", "Swap"};
  
  Forest mu_forest(n_tree_mu, 1, n);
  
  Forest growth_forest(n_tree_growth, 1, n);
  
  Forest treat_forest(n_tree_treat, 1, n);
  
  //for holding covariate selections
  arma::cube covariates_treat(p_treat, n_tree_treat, n_iter);
  
  //for holding covariate selections
  arma::cube covariates_growth(p_growth, n_tree_growth, n_iter);
  
  //initialise missing z values
  for(int i = 0; i<n; i++)
  {
    if(z_missing[i])
    {
      if(p_scores[i]>=0.5)
      {
        Z_treat[i]=1;
      }
      else
      {
        Z_treat[i]=0;
      }
    }
  }
  
  //Matrix to store if changed or not changed
  LogicalMatrix z_changed(sum(z_missing), num_keep);
  
  for(int iter = 0; iter < n_iter; iter++)
  {
    for(int tree_num = 0; tree_num < n_tree_mu; tree_num++)
    {
      NumericVector y_resid = y_scaled-rowSumsWithoutColumn(tree_preds_mu, tree_num)-Z_growth*rowSumsWithoutColumn(tree_preds_growth, -1)-Z_treat*rowSumsWithoutColumn(tree_preds_treat, -1);
      
      String choice = sample(choices, 1)[0];
      
      Tree proposal_tree = Tree(mu_forest.tree_vector[tree_num]);
      
      if(choice == "Grow")
      {
        proposal_tree.grow(X_mu, p_mu, min_nodesize);
      }
      
      
      if(choice == "Prune")
      {
        proposal_tree.prune();
      }
      
      
      if(choice == "Change")
      {
        proposal_tree.change(X_mu, p_mu);
        proposal_tree.change_update(X_mu);
      }
      
      
      if(choice == "Swap")
      {
        proposal_tree.swap();
        proposal_tree.change_update(X_mu);
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
      
      for(int i=0; i<n; i++)
      {
        tree_preds_mu(i, tree_num) = tree_preds_from_iter_mu[i];
      }
    
    }
    
    for(int tree_num = 0; tree_num < n_tree_growth; tree_num++)
    {
      NumericVector y_resid = y_scaled-rowSumsWithoutColumn(tree_preds_mu, -1)-Z_growth*rowSumsWithoutColumn(tree_preds_growth, tree_num)-Z_treat*rowSumsWithoutColumn(tree_preds_treat, -1);
      
      String choice = sample(choices, 1)[0];
      
      Tree proposal_tree = Tree(growth_forest.tree_vector[tree_num]);
      
      if(choice == "Grow")
      {
        proposal_tree.grow(X_growth, p_growth, min_nodesize);
      }
      
      
      if(choice == "Prune")
      {
        proposal_tree.prune();
      }
      
      
      if(choice == "Change")
      {
        proposal_tree.change(X_growth, p_growth);
        proposal_tree.change_update(X_growth);
      }
      
      
      if(choice == "Swap")
      {
        proposal_tree.swap();
        proposal_tree.change_update(X_growth);
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
      
      for(int i=0; i<n; i++)
      {
        tree_preds_growth(i, tree_num) = tree_preds_from_iter_growth[i];
      }
      
      for(int t = 0; t<growth_forest.tree_vector[tree_num].node_vector.size(); t++)
      {
        if(growth_forest.tree_vector[tree_num].node_vector[t].in_use & !growth_forest.tree_vector[tree_num].node_vector[t].is_terminal)
        {
          int chosen_var = growth_forest.tree_vector[tree_num].node_vector[t].variable;
          
          covariates_growth(chosen_var, tree_num, iter) += 1;
        }
      }
      
    }
    
    
    for(int tree_num = 0; tree_num < n_tree_treat; tree_num++)
    {
      NumericVector y_resid = y_scaled-rowSumsWithoutColumn(tree_preds_mu, -1)-Z_growth*rowSumsWithoutColumn(tree_preds_growth, -1)-Z_treat*rowSumsWithoutColumn(tree_preds_treat, tree_num);
      
      String choice = sample(choices, 1)[0];
      
      Tree proposal_tree = Tree(treat_forest.tree_vector[tree_num]);
      
      if(choice == "Grow")
      {
        proposal_tree.grow(X_treat, p_treat, min_nodesize);
      }
      
      
      if(choice == "Prune")
      {
        proposal_tree.prune();
      }
      
      
      if(choice == "Change")
      {
        proposal_tree.change(X_treat, p_treat);
        proposal_tree.change_update(X_treat);
      }
      
      
      if(choice == "Swap")
      {
        proposal_tree.swap();
        proposal_tree.change_update(X_treat);
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
      
      for(int i=0; i<n; i++)
      {
        tree_preds_treat(i, tree_num) = tree_preds_from_iter_treat[i];
      }
      
      for(int t = 0; t<treat_forest.tree_vector[tree_num].node_vector.size(); t++)
      {
        if(treat_forest.tree_vector[tree_num].node_vector[t].in_use & !treat_forest.tree_vector[tree_num].node_vector[t].is_terminal)
        {
          int chosen_var = treat_forest.tree_vector[tree_num].node_vector[t].variable;
          
          covariates_treat(chosen_var, tree_num, iter) += 1;
        }
      }
      
    }
      
    Rcpp::Rcout << "Total of " << iter+1 << " of " << n_iter << " iterations completed! " << "(" << (float)(iter+1)/(float)n_iter*100 << "%)             " << "\r";
    Rcpp::Rcout.flush();
    
    NumericVector iter_preds_mu = rowSumsWithoutColumn(tree_preds_mu, -1);
    NumericVector iter_preds_growth = rowSumsWithoutColumn(tree_preds_growth, -1);
    NumericVector iter_preds_treat = rowSumsWithoutColumn(tree_preds_treat, -1);
    
    
    if(iter>=(n_burn-1) & (iter-n_burn)%keep_every==0)
    {
    for(int i=0; i<n; i++)
    {
      preds_mat_mu(i, num_kept) = iter_preds_mu[i];
      preds_mat_growth(i, num_kept) = iter_preds_growth[i];
      preds_mat_treat(i, num_kept) = iter_preds_treat[i];
    }
    }
    
    tau=sample_tau(n, nu, y_scaled, iter_preds_mu+Z_growth*iter_preds_growth+Z_treat*iter_preds_treat, lambda);
    
    if(iter>=(n_burn-1) & (iter-n_burn)%keep_every==0)
    {
    taus[num_kept] = tau;
    }
    
    
    if(iter>=(n_burn-1) & (iter-n_burn)%keep_every==0)
    {
    int my_counter = 0;
    //impute missing z values again
    for(int i = 0; i<n; i++)
    {
      if(z_missing[i])
      {
        double new_z_prob = prob_z(p_scores[i], tau, y_scaled[i], iter_preds_mu[i], iter_preds_growth[i], iter_preds_treat[i]);
        
        double runif_value = R::runif(0, 1);
        
        if(new_z_prob>=runif_value)
        {
          Z_treat[i]=1;
          
          z_changed(my_counter, num_kept) = true;
        }
        else
        {
          Z_treat[i]=0;
          
          z_changed(my_counter, num_kept) = false;
        }
        my_counter++;
      }
    }
    num_kept = num_kept + 1;
    }
  }
  
  Rcpp::Rcout << "";
  
  
  return List::create(
    Named("predictions_mu") = (preds_mat_mu*y_sd)+y_mean,
    Named("predictions_growth") = preds_mat_growth*y_sd,
    Named("predictions_treat") = preds_mat_treat*y_sd,
    Named("taus") = taus/(pow(y_sd, 2)),
    Named("sigmas") = y_sd/sqrt(taus),
    Named("changes") = z_changed,
    Named("covariates_growth") = covariates_growth,
    Named("covariates_treat") = covariates_treat
  );
}
 

