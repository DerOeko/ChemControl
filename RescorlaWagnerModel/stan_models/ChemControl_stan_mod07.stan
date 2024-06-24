data {
    int<lower=1> nTrials;
    int<lower=1> nBlocks;
    int<lower=1> nSubjects;
    int<lower=1, upper=4> state[nSubjects, nBlocks, nTrials];
    int<lower=0, upper=1> action[nSubjects, nBlocks, nTrials];
    int<lower=-1, upper=1> feedback[nSubjects, nBlocks, nTrials];
}

transformed data {
   vector[4] initQ; 
   initQ = [0.5, -0.5, 0.5, -0.5]';
}

parameters {
   vector[7] mu_pr;
   vector<lower=0>[7] sigma;
   matrix[nSubjects, nBlocks] ep_pr;
   matrix[nSubjects, nBlocks] rho_pr;
   matrix[nSubjects, nBlocks] goBias_pr;
   matrix[nSubjects, nBlocks] omegaInit_pr;
   matrix[nSubjects, nBlocks] alphaO_pr;
   matrix[nSubjects, nBlocks] betaO_pr;
   matrix[nSubjects, nBlocks] thresO_pr;
}

transformed parameters {
    matrix<lower=0, upper=1>[nSubjects, nBlocks] ep;
    matrix<lower=0>[nSubjects, nBlocks] rho;
    matrix[nSubjects, nBlocks] goBias;
    matrix<lower=0, upper=1>[nSubjects, nBlocks] omegaInit;
    matrix<lower=0, upper=1>[nSubjects, nBlocks] alphaO;
    matrix<lower=0>[nSubjects, nBlocks] betaO;
    matrix<lower=-1, upper=1>[nSubjects, nBlocks] thresO;

    // Generate epsilon values based on priors.
    for (s in 1:nSubjects) {
        for (b in 1:nBlocks){
            ep[s, b]  = Phi_approx(mu_pr[1] + sigma[1] * ep_pr[s, b]);
            omegaInit[s, b] = Phi_approx(mu_pr[4] + sigma[4] * omegaInit_pr[s, b]);
            alphaO[s, b] = Phi_approx(mu_pr[5] + sigma[5] * alphaO_pr[s,b]);
            thresO[s, b] = tanh(0.5*(mu_pr[7] + sigma[7] * thresO_pr[s, b]));
        }
    }
        
    // Confine rho to a strictly positive range (exp() does that)
    rho = exp(mu_pr[2] + sigma[2] * rho_pr);
    goBias = mu_pr[3] + sigma[3] * goBias_pr;
    betaO = exp(mu_pr[5] + sigma[5] * betaO_pr);
 }

model {
    mu_pr[1] ~ normal(0, 1.0);
    mu_pr[2] ~ normal(0, 1.0);
    mu_pr[3] ~ normal(0, 10.0);
    mu_pr[4] ~ normal(0, 1.0);
    mu_pr[5] ~ normal(0, 1.0);
    mu_pr[6] ~ normal(0, 1.0);
    mu_pr[7] ~ normal(0, 1.0);

    sigma[1] ~ normal(0, 0.2);
    sigma[2] ~ normal(0, 0.2);
    sigma[3] ~ cauchy(0, 1.0);
    sigma[4] ~ normal(0, 0.2);
    sigma[5] ~ normal(0, 0.2);
    sigma[6] ~ normal(0, 0.2);
    sigma[7] ~ normal(0, 0.2);


    for (s in 1:nSubjects) {
        ep_pr[s] ~ normal(0, 1.0);
        rho_pr[s] ~ normal(0, 1.0);
        goBias_pr[s] ~ normal(0, 1.0);
        omegaInit_pr[s] ~ normal(0, 1.0);
        alphaO_pr[s] ~ normal(0, 1.0);
        betaO_pr[s] ~ normal(0, 1.0);
        thresO_pr[s] ~ normal(0, 1.0);

        for (b in 1:nBlocks){
            vector[4] w_g;  // action weight for go
            vector[4] w_ng; // action weight for nogo
            vector[4] q_g;  // Q value for go
            vector[4] q_ng; // Q value for nogo
            vector[4] pGo;   // prob of go (press)
            vector[4] sv;
            real Omega;
            real omega;
            real q_pe;
            real v_pe;

            w_g  = initQ*rho[s, b];
            w_ng = initQ*rho[s, b];
            q_g  = initQ*rho[s, b];
            q_ng = initQ*rho[s, b];
            sv = initQ;
            Omega = 0;
            omega = omegaInit[s, b];

            for (t in 1:nTrials){
                w_g[state[s, b, t]] = (1-omega) * q_g[state[s, b, t]] + goBias[s, b] + omega * sv[state[s, b, t]];
                w_ng[state[s, b, t]] = (1-omega) * q_ng[state[s, b, t]];

                action[s, b, t] ~ bernoulli_logit(w_g[state[s, b, t]] - w_ng[state[s, b, t]]);

                v_pe = rho[s, b] * feedback[s, b, t] - sv[state[s, b, t]];
                sv[state[s,b,t]] += ep[s, b] * (rho[s, b] * feedback[s, b, t] - sv[state[s, b, t]]);

                if (action[s, b, t]) { // update go value
                    q_pe = rho[s, b] * q_g[state[s,b,t]];
                    q_g[state[s, b, t]] += ep[s, b] * (rho[s, b] * feedback[s, b, t] - q_g[state[s, b, t]]);
                } else { // update no-go value
                    q_pe = rho[s, b] * q_ng[state[s,b,t]];
                    q_ng[state[s, b, t]] += ep[s, b] * (rho[s, b] * feedback[s, b, t] - q_ng[state[s, b, t]]);
                } // end of if clause

                Omega = Omega + alphaO[s, b]*(q_pe - v_pe - Omega);
                omega = 1/(1+exp(-betaO[s, b]*(Omega-thresO[s, b])));
                
            } // end of t loop
        } // end of b loop
    } // end of s loop
} // end of model block