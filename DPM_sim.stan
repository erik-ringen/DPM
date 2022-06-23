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

} // end of functions block


data{
  int N; // number of tips
  int J; // number of response traits
  int resp_type[J]; // type of trait (1 = gaussian, 2 = binomial)
  int N_seg; // total number of segments in the tree
  int node_seq[N_seg]; // index of tree nodes
  int parent[N_seg]; // index of the parent node of each descendent
  vector[N_seg] ts; // time since parent
  vector[N_seg] tip; // indicator of whether a given segment ends in a tip
}

parameters{
vector[J] b; // SDE intercepts
vector[J] eta_anc; // ancestral states
matrix[N_seg - 1,J] z_drift; // stochastic drift, unscaled and uncorrelated

real a_obs[J]; // intercept for observed variables
}

transformed parameters{
matrix[N_seg,J] eta;
matrix[J,J] Q; // "drift" matrix
matrix[J,J] I; // identity matrix
matrix[J,J] A;  // "selection" matrix
matrix[N_seg,J] drift_tips; // terminal drift parameters, saved here to use in likelihood for Gaussian outcomes
matrix[N_seg,J] sigma_tips; // terminal drift parameters, saved here to use in likelihood for Gaussian outcomes

// selection matrix, col 1 is var1, col 2 is var2
A[1,1] = -0.2; // autoregressive effect of var1 on itself
A[2,2] = -0.2; // autoregressive effect of var2 on itself

A[2,1] = 4;
A[1,2] = 0;

// drift matrix
Q[1,1] = 2;
Q[2,2] = 2;
Q[1,2] = 0;
Q[2,1] = 0;

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
}

generated quantities{
  real y[N,2];

  for (j in 1:J) {
  
    if (resp_type == 1) {
      for (i in 1:N) y[i,j] = normal_rng( eta[i,j], sigma_tips[i,j] );
    }
    
    if (resp_type == 2) {
      for (i in 1:N) y[i,j] = bernoulli_logit_rng( eta[i,j] + drift_tips[i,j] );
    }
}
}

