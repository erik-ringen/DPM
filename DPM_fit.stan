functions { 
matrix kronecker_prod(matrix A, matrix B) { // returns the Kronecker Product
  matrix[rows(A) * rows(B), cols(A) * cols(B)] C;
  int m;
  int n;
  int p;
  int q;
  m = rows(A);
  n = cols(A);
  p = rows(B);
  q = cols(B);
  for (i in 1:m)
    for (j in 1:n)
      for (k in 1:p)
        for (l in 1:q)
      C[p*(i-1)+k,q*(j-1)+l] = A[i,j]*B[k,l];
  return C;
}

matrix A_dt(matrix A, real t) {  // expected auto and cross effects over a discrete time t
 return( matrix_exp(A * t) );
}

matrix A_sharp(matrix A) {
  matrix[rows(A) * rows(A), cols(A) * cols(A)] A_temp;
  matrix[rows(A),cols(A)] I; // identity matrix
  
  I = diag_matrix(rep_vector(1,rows(A)));

  A_temp = kronecker_prod(A,I) + kronecker_prod(I,A);
  return(A_temp);
}

matrix cov_drift(matrix A, matrix Q, real ts) {
  matrix[rows(A) * rows(A), cols(A) * cols(A)] A_sharp_temp;
  matrix[rows(A) * rows(A), cols(A) * cols(A)] I; // identity matrix
  vector[rows(Q)*cols(Q)] row_Q;
  vector[rows(A)*cols(A)] irow_vec;
  matrix[rows(A),cols(A)] irow_mat;
  
  I = diag_matrix(rep_vector(1,rows(A_sharp_temp)));
  
  A_sharp_temp = A_sharp(A);
  
  // row operation takes elements of a matrix rowwise and puts them into a column vector
  for (i in 1:rows(Q))
  for (j in 1:cols(Q)) {
    row_Q[i + (j-1)*rows(Q)] = Q[j,i];
  }
  
  irow_vec = inverse(A_sharp_temp) * (matrix_exp(A_sharp_temp * ts) - I) * row_Q;
  
  // irow takes elements of a column vector and puts them in a matrix rowwise
  {
  int row_size = rows(A);
  int row_ticker = 1;
  int col_ticker = 0;
  
  for (i in 1:num_elements(irow_vec)) {
    col_ticker += 1;
    if (col_ticker > row_size) {
      row_ticker += 1;
      col_ticker = 1;
    }
    irow_mat[row_ticker,col_ticker] = irow_vec[i];
  }
  }
  return(irow_mat);
}

// Function to convert from type real to type interger,when appropriate
 int bin_search(real x, int min_val, int max_val){
    // This assumes that min_val >= 0 is the minimum integer in range, 
    //  max_val > min_val,
    // and that x has already been rounded. 
    //  It should find the integer equivalent to x.
    int range = (max_val - min_val+1) %/% 2; // We add 1 to make sure that truncation doesn't exclude a number
    int mid_pt = min_val + range;
    int out;
    while(range > 0) {
      if(x == mid_pt){
        out = mid_pt;
        range = 0;
      } else {
        // figure out if range == 1
        range =  (range+1) %/% 2; 
        mid_pt = x > mid_pt ? mid_pt + range: mid_pt - range; 
        }
    }
    return out;
  }
} // end of functions block


data{
  int N_tips; // number of tips
  int J; // number of response traits
  array[J] int resp_type; // type of trait (1 = gaussian, 2 = binomial)
  int N_seg; // total number of segments in the tree
  array[N_seg] int node_seq; // index of tree nodes
  array[N_seg] int parent; // index of the parent node of each descendent
  vector[N_seg] ts; // time since parent
  vector[N_seg] tip; // indicator of whether a given segment ends in a tip
  array[N_tips, J] real X;
}

parameters{
vector<upper=0>[J] A_diag; // auto-regressive terms of A
vector[J*J - J] A_offdiag; // 
vector<lower=0>[J] Q_diag; // self drift terms
vector[J] b; // intercepts ("bias")
vector[J] eta_anc; // ancestral states
matrix[N_seg - 1,J] z_drift; // stochastic drift, unscaled and uncorrelated
}

transformed parameters{
matrix[N_seg,J] eta;
matrix[J,J] Q; // "drift" matrix
matrix[J,J] I; // identity matrix
matrix[J,J] A;  // "selection" matrix
matrix[N_seg,J] drift_tips; // terminal drift parameters, saved here to use in likelihood for Gaussian outcomes
matrix[N_seg,J] sigma_tips; // terminal drift parameters, saved here to use in likelihood for Gaussian outcomes
vector[J] Q_offdiag = rep_vector(0.0, J);

// Fill A matrix //////////
{
int ticker = 1;
// fill upper tri of matrix
for (i in 1:(J-1))
  for (j in (i+1):J) {
    A[i,j] = A_offdiag[ticker];
    ticker += 1;
  }
// fill lower tri of matrix
  for (i in 1:(J-1))
  for (j in (i+1):J) {
    A[j,i] = A_offdiag[ticker];
    ticker += 1;
  }
// fill diag of matrix
  for (j in 1:J) A[j,j] = A_diag[j];
}

// Fill Q matrix /////
{
int ticker = 1;
for (i in 1:(J-1))
for (j in (i+1):J) {
  Q[i,j] = Q_offdiag[ticker];
  Q[j,i] = Q_offdiag[ticker]; // symmetry of covariance
  ticker += 1;
}
for (j in 1:J) Q[j,j] = Q_diag[j];
}

// identity matrix
I = diag_matrix(rep_vector(1,J));

// setting ancestral states, and placeholders
for (j in 1:J) {
eta[node_seq[1],j] = eta_anc[j]; // ancestral state
drift_tips[node_seq[1],j] = -99; // 
sigma_tips[node_seq[1],j] = -99; // 
}

  for (i in 2:N_seg) {
  matrix[J,J] A_delta; // amount of deterministic change ("selection")
  matrix[J,J] VCV; // variance-covariance matrix of stochastic change ("drift")
  vector[J] drift_seg; // accumulated drift over the segment
  
  A_delta = A_dt(A, ts[i]);
  VCV = cov_drift(A, Q, ts[i]);
  
  // No drift on the interaction, bc its simply a product of RI and TPC
  drift_seg = cholesky_decompose(VCV) * to_vector( z_drift[i-1,] );
  
  // if not a tip, add the drift parameter
  if (tip[i] == 0) {
  eta[node_seq[i],] = to_row_vector(
    A_delta * to_vector(eta[parent[i],]) + (inverse(A) * (A_delta - I) * b) + drift_seg
  );

  drift_tips[node_seq[i],] = to_row_vector(rep_vector(-99, J));
  sigma_tips[node_seq[i],] = to_row_vector(rep_vector(-99, J));

  }
  // if is a tip, omit, we'll deal with it in the model block
  else {
    eta[node_seq[i],] = to_row_vector(
    A_delta * to_vector(eta[parent[i],]) + (inverse(A) * (A_delta - I) * b)
  );

  drift_tips[node_seq[i],] = to_row_vector(drift_seg);
  sigma_tips[node_seq[i],] = to_row_vector(sqrt(diagonal(VCV)));
  }

  }

} // end transformed parameters block

model{
  b ~ std_normal();
  eta_anc ~ std_normal();
  to_vector(z_drift) ~ std_normal();
  A_offdiag ~ std_normal();
  A_diag ~ std_normal();
  Q_diag ~ std_normal();
  
  for (j in 1:J) {
      if (resp_type[j] == 1) {
        for (i in 1:N_tips) (X[i,j] - eta[i,j]) ~ normal(0, sigma_tips[i,j]);
    }
      if (resp_type[j] == 2) {
        for (i in 1:N_tips) bin_search(X[i,j], 0, 1) ~ bernoulli_logit(eta[i,j] + drift_tips[i,j]);
    }
  }
}
