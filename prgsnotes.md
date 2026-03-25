## Model specification

Let $X_{ij}$ be the number of times agent $i$ beats agent $j$ out of the $n_{ij}$ total dyadic interactions between the agents. We model $X_{ij}$ to be binomially distributed with $p_{ij}$ being the probability that $i$ beats $j$ in a single interaction:

$X_{ij} \sim \text{Bin}(n_{ij},p_{ij})$.

To express intransitivity in the model, we impose a certain bilinear form on the logit of the binomial probability

$\text{logit }p_{ij}=a(z_i-z_j) + b(y_i-y_j) + c(z_iy_j-z_jy_i)$,

where $y_i$ is an observed covariate describing an aspect of the dominance of agent $i$, while $z_i$ is a latent variable, representing some unobserved characteristic of agent $i$. 

## Exploration

- Investigate which combinations of interaction count intensity parameter $\lambda$ and agent count $s$ makes for unstable inference. How few dyadic interactions is possible to still be able to get some useful inference?
- Predicitve check between intransitive model and transitive model with $c=0$. 
- Investigate what priors are good for model parameters. Sensitivity analysis. Informative vs. domain specific priors. 
- Fit the model on real data. 
- Count number of triadic intransitive cycles, and compare to expected number of cycles from model binomial probabilities.
- Expand to time series dominance. Evolution of model parameters over time. 
- Comparing ranking from model to traditional ranking methods such as Davis's score. 
- Identifiability, posterior geometry
- Når begge (y,z) er latente, kan forankre den ene y i null f.eks. for å forhindre ikke identifiserbarhet gjennom rotasjonssymmetri. 
- ARIMA-modeller: Hvis man får et svar som er utenfor domenet, sett til en verdi som er innenfor domenet ved symmetrier. 

- Er modellen bedre enn standard Bradley Terry modell
- gjøre modellen en naturlig utvidelse av R bradely terry. 
- bevise hvordan symmetritransformasjoner fungerer. 
- Likelihood ratio test, intransitivitet mot ikke transitivitet. 
- Store datasett der man kan tenke seg at det er transitivitet.
- gjøre om modell på matrisevektorform. 

- forskjell mellom mye data og lite data (antall diader)

- Dominanshierarkier over tid: er de stabile? 

- Skrive om David's score
