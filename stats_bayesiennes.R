###########################################################################
####### CODE DE L'UTILISATION DE L'ALGORITHME SUR DONNEES SIMULEES ########
###########################################################################



# Chargement des libraries ------------------------------------------------

library("statmod")
library("sda")
library("GIGrvg")
library("MCMCpack")
library("mvtnorm")


# Création d’une fonction pour l'algorithme -------------------------------

# La fonction prend en inputs:
# n : la dimension du paramètre à estimer
# qn: le niveau de parcimonie
# N_iter: le nombre d'itérations pour calculer obtenir la distribution a posteriori
# N_replications : le nombre de répétitions de l'algorithme entier (simulation de y comprise)

# L'output est le suivant:
# Affiche la somme des erreurs quadratiques entre médianes et vrais paramètres
# Renvoie med : le vecteur des médianes a posteriori des n paramètres

gibbs <- function(n,qn,A,N_iter,N_replications){
  # On cree une matrice de stockage de resultats
  results = matrix(data = NA, nrow = N_replications, ncol = n)
  #chaque ligne correspond à une répétition de la simulation et chaque colonne j à un theta_j
  
  # On cree une matrice de stockage des médianes
  med = rep(1,n)
  
  # On crée une boucle pour le nombre de répétitions de l'algorithme
  for (n_replications in 1:N_replications){
    theta = rep(1, n)
    psi = rep(1, n)
    phi = rep(1, n)
    a = rep(1/2,n)
    tau = 1
    t = rep(1, n)
    
    # On simule l'ensemble des priors
    psi = rexp(n = n, rate = 1/2)
    phi = rdirichlet(n = 1, alpha = a)
    tau = rgamma(n = 1, shape = n*a, scale = 1/2) #scale ou rate?
    y = rmvnorm(1, mean=c(rep(A, n*qn),rep(0, n-n*qn)), sigma=diag(n))
    a=1/2
    
    mean = 0
    # On calcule les posteriors en itérant l'algorithme N_iter fois
    for(n_iterations in 1:N_iter){
      # Simulation des theta
      for (j in 1:n){
        theta[j] = rnorm(1, mean= y[j]*((1+1/(psi[j]*phi[j]^ 2*tau^2))^(-1)), sd= sqrt((1+1/(psi[j]*phi[j]^2*tau^2))^(-1)))
      }
      # Simulation des psi
      for (j in 1:n){
        psi[j] = (rinvgauss(1, mean = phi[j]*tau/abs(theta[j]), shape = 1))
      }
      # Simulation de tau
      tau = rgig(n = 1, lambda = n*a-n, psi = 1, chi = 2*sum(abs(theta)/phi))
      # Simulation des T
      for(j in 1:n){
        t[j] = rgig(n = 1,  lambda = a-1, psi = 1, chi = 2*abs(theta[j]))
      }
      # Simulation des phi
      for(j in 1:n){
        phi[j] = t[j]/sum(t)
      }
      burnin = N_iter/2
      # On ne conserve que les theta après la période de burn-in
      mean = mean + (n_iterations > burnin)*theta/(N_iter-burnin)
    }
    
    # On sauvegarde avant de passer a la répétition suivante de tout l'algorithme
    for(j in 1:n){
      results[n_replications,j] = mean[j]
    }
    # On met en place un compeur pour suivre l'avancement de l'algorithme
    print(n_replications)
  }
  
  # Une fois les répétitions faites, on calcule ensuite la médiane 
  # de chaque colonne, càd de chaque theta_j
  for(j in 1:n){
    med[j]=median(results[,j])
  }
  print("Resultat de la somme de l'erreur quadratique")
  print(sum((med-c(rep(A, n*qn),rep(0, n-n*qn)))^2))
  return(med)
}


# Mise en place des differents scenarii -----------------------------------

# Exemple de lancement d'un scenario
med = gibbs(n = 100, qn = 0.05, A = 7, N_iter = 3000, N_replications = 100)

# On recupere la mediane et on l'affiche
print(med)

# On affiche les différentes médianes
plot(x = 1:100, y = med, xlab = "Indice du paramètre", ylab = "Médiane", 
     main = "Médiane a posteriori avec prior Di(1/2)", type = "p")






###################### Application au cancer de la prostate ##############



###### Ouvertude du jeu de données ##############
data(singh2002)
hist(singh2002$x)
singh2002$y
df = t(singh2002$x)
colnames(df) <- c(rep(x = "healthy", 50),rep(x = "cancer", 52)) 
#df est une base de données contenant 102 individus (50 sains, 52 avec le cancer) et 6033 genes pour chaque individu
sains = as.data.frame(df[,1:50])
cancer = as.data.frame(df[,51:102])


# t test ------------------------------------------------------------------
#On va réaliser un two sample t-test dans un premier temps
#test pour une ligne 
test = t.test(sains[1,],cancer[1,], var.equal=TRUE, paired=FALSE)

#test multidimensionnel
test=matrix(NA, ncol = 1, nrow = 6033)

for (i in 1:6033){
  test[i,1] = t.test(sains[i,],cancer[i,], var.equal=TRUE, paired=FALSE)$statistic
}

#Plot
hist(test[,1], freq = F, breaks = 60, ylim = c(0, 0.4))
curve(dnorm(x, 0 , 1), col="red", add=T)

#On passe à la z-stat
z = qnorm(pt(test, df = 100))
#graph de la stat
hist(z[,1], freq = F, breaks = 40, ylim = c(0, 0.4), main = "histogramme des z-values", xlab = "")
curve(dnorm(x, 0 , 1), col="red", add=T)

head(z[,1])

# Algorithme --------------------------------------------------------------
#Algorithme (partie 2.4)
#Initialisation
n=6033
N_iter = 10000
theta = rep(1, n)
psi = rep(1, n)
phi = rep(1, n)
y = df[1,]
tau = 1
#a=runif(n, 1/1000, 1/2)
#a = rep(1/20,n)
a = rep(runif(1, 1/6000, 1/2),n)

m1=matrix(NA, nrow = N_iter, ncol = n)
t = vector(mode = "numeric", length = n)


# On simule l'ensemble des priors
psi = rexp(n = n, rate = 1/2)
phi = rdirichlet(n = 1, alpha = a)
tau = rgamma(n = 1, shape = n*a, scale = 1/2) #scale ou rate?
for (j in 1:n){
  theta[j] = rnorm(1, mean = 0, sd = sqrt(psi[j]*phi[j]^2*tau^2))
}

for(i in 1:N_iter){
  a = rep(runif(1, 1/6000, 1/2),n)
  # Simulation des theta
  for (j in 1:n){
    theta[j] = rnorm(1, mean= z[j]*((1+1/(psi[j]*phi[j]^2*tau^2))^(-1)), sd= sqrt((1+1/(psi[j]*phi[j]^2*tau^2))^(-1)))
  }
  # Simulation des psi
  for (j in 1:n){
    psi[j] = rinvgauss(1, mean = phi[j]*tau/abs(theta[j]), shape = 1)
  }
  # Simulation de tau
  #tau = rgig(n = 1, lambda = lambda-n, psi = 1, chi = 2*sum(abs(theta)/phi))
  tau = rgig(n = 1, lambda = n*a-n, psi = 1, chi = 2*sum(abs(theta)/phi))
  # Simulation des T
  for(j in 1:n){
    t[j] = rgig(n = 1,  lambda = a-1, psi = 1, chi = 2*abs(theta[j]))
  }
  # Simulation des phi
  for(j in 1:n){
    phi[j] = t[j]/sum(t)
  }
  m1[i,]=theta
}


# clustering --------------------------------------------------------------

#initialisation
burn_in = 5000
cluster = matrix(NA, nrow = N_iter-burn_in, ncol = 6033)
cluster_center = matrix(NA, nrow = N_iter-burn_in, ncol = 2)
cluster_size = matrix(NA, nrow = 2, ncol = N_iter-burn_in)
min_cluster_size = rep(1, n)

#boucle
for(i in (burn_in+1):N_iter){
  k=kmeans(abs(m1[i,]), centers = 2)
  cluster[i-burn_in,]=k$cluster
  cluster_size[,i-burn_in] = k$size
  cluster_center[i-burn_in,] = k$center
}

hist(cluster_size, breaks = 1000, freq = F, xlab = "Taille des clusters", main = "Histogramme de la taille des clusters")

#for (j in 1:N_iter-burn_in){
#  min_cluster_size[j] = 6033 - max(cluster_size[1,j],cluster_size[2,j], cluster_size[3,j])
#}


for (j in 1:N_iter-burn_in){
  min_cluster_size[j] = min(cluster_size[1,j],cluster_size[2,j])
}
summary(min_cluster_size)
min_cluster_size

#On cherche le mode 
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
getmode(min_cluster_size)

### Extraction des gênes les plus signigicatifs 

#Avec les gènes dont le niveau diffère le plus
gene = matrix(NA, nrow = 6033, ncol = 1)
for (i in 1:6033){
  gene[i,1] = median(abs(m1[5001:10000,i]))
}
gene

#Convergence du gibbs pour un theta donné : le cas du gène 601
plot(cumsum(abs(m1[1:10000, 601]))/1:10000, type = "l", ylab = "", xlab = "itérations", main = "Convergence de l'estimateur du gène 601")
abline(v = 5000, col = "red")