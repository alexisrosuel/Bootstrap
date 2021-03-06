#' ---
#' title: "Projet Applications du bootstrap et autres techniques de r��chantillonage"
#' author: "Pascal Jauffret et Alexis Rosuel"
#' date: " 28 Avril 2017"
#' output: pdf_document
#' ---

#' Au cours de ce projet, nous avons d�cid� de nous interesser � un jeu de donn�es fournit par le Texas Department of Justice, disponible ici :
#' http://www.tdcj.state.tx.us/death_row/dr_executed_offenders.html. Ce dataset rassemble quelques informations basiques concernant
#' tous les condamn�s � la peine de mort et �xecut�s depuis le 7 d�cembre 1982, jusqu'au 14 mars dernier. En particulier on connait
#' leur �ge, la date de l'ex�cution et leur ethnie.


#' 
#' Description du jeu de donn�es
#' ===

#' Commen�ons par importer le jeu de donn�es
setwd('C:\\Users\\Alexis\\Google Drive\\Documents\\ENSAE\\Semestre 2\\Applications du bootstrap\\examen')
data = read.csv('executed.csv', sep=';')

#' Affichons les 6 premi�res lignes
head(data)

#' Par la suite, on va r�pondre aux questions suivante :  
#' - la distribution de l'�ge des condamn�s est-elle issue d'une loi normale ?
#' - quelle est leur moyenne d'�ge ?  
#' - le jour de leur ex�cution est-elle uniforme ? (ou bien il y en a plus le mardi que le dimanche par exemple ?)


#' 
#' Mod�lisation et justification th�orique des hypoth�ses n�cessaires
#' ===
#' 
#' Afin d'appliquer les m�thodes de r��chantillonage vues en cours, nous avons besoin de fixer un cadre th�orique � notre �tude.
#' Notons $X_i \in R^k$ chaque individu. Ici, on a pas de raison de penser qu'il puisse y avoir de lien entre les individus, ni
#' qu'ils ne soient pas issus de la m�me distribution, qu'on notera $P$.  
#' 
#' Pour r�sumer, on a donc ceci :
#' $\mathcal{X}_n = (X_1, ..., X_n) \sim P$ iid.
#' 
#' L'�ge est une variable al�atoire born�e, tout comme le jour de l'�xecution. Ainsi, tous les moments des 
#' composantes des individus $X_i$ sont bien d�finis. L'in�galit� de Hoeffding pour les U-statistiques est donc aussi applicable, 
#' ce qui nous permet d'appliquer le subsampling.
#' 
#' $$E[\|X_n\|] < +\infty \text{ et } E[\|X_n\|^2] < +\infty$$
#' 
#' Par ailleurs, la condition (ii) du th�or�me de Lindeberg bien v�rifi�e. 
#' 
#' Enfin, tous nos estimateur (moyenne, min, max, etc.) seront bien invariants par permutation.



#' 
#' Simulations
#' ===
#' 
#' On charge la librairie 'boot' qui permet d'appliquer le bootstrap de mani�re efficace dans R
library(boot)

#' Dans tout ce qui suit, on utilisera B = 10000 simulations
B = 10000

#' On dispose de 542 observations dans notre dataset :
n = length(data$Age)
n

#' Et on sous �chantillonera des sous ensembles de taille b = 25 (on prend en fait $b \approx \sqrt{n}$)
b = 25

#' Concernant les tests, on se placera toujours au niveau $\alpha=0.05$
alpha = 0.05

#' Enfin, et afin de pouvoir reproduire les r�sultats, on fixe :
set.seed(123)


#' 
#' Application
#' ===
#' 
#' Etude de la variable Age
#' ---
#' 
#' Dans cette section, nous allons nous concentrer sur la variable �ge. Les questions auxquelles nous allons essayer de r�pondre
#' sont celles-ci:  
#' 1) Est ce que la distribution des �ges est normale ?  
#' 1) Quel est un intervalle de confiance de la moyenne des condamn�s au niveau 95% ? (bootstrap vs subsampling vs test de Fisher)    
#' 2) L'�ge moyen est-il de 38 ans ? (bootstrap sous contrainte et par symm�trisation)    
#' 
#' Commen�ons par le test de Kolmogorov-Smirnov en utilisant la m�thode du bootstrap param�trique. A premi�re vue, on s'attend �
#' r�pondre que la distribution r�elle n'est pas issue d'une loi normale :
hist(data$Age, nclass=25, freq=FALSE)
curve(dnorm(x,mean(data$Age),sd(data$Age)), add=TRUE, col='blue')

#' Plus pr�cisiement, on note $\hat{F}_n$ la fonction de r�partition empirique associ�e. 
#' On note aussi $\mathcal{P}_0 = \{\mathcal{N}(\mu, \sigma^2), \mu \in R, \sigma^2 \in R^+\}$ la famille de lois dont on souhaite 
#' tester si nos donn�es sont issues.
#' 
#' La statistique de test naturelle est $T_n(\mathcal{X}_n) =  \sqrt{n} \|\hat{F}_n - F_0 \|_\infty$. Enfin, on note la loi de
#' $T_n(\mathcal{X}_n)$ : $L_n(P_0,P_0) = \mathcal{L}(\sqrt{n} \|\hat{F}_n - F_0 \|_\infty)$
#' 
#' On calcule $\hat{\theta}_n$ l'estimateur de $\theta$ : 
mu = mean(data$Age)
sigma2 = sd(data$Age)

#' On calcule la statistique $T_n$ :
norme_infinie <- function(Xn){
  Xn.sorted = sort(Xn)
  resultat = c()
  n = length(Xn)
  for(i in 1:n){
    resultat = c(resultat ,max(i/n - pnorm(Xn.sorted[i], mean=mu, sd=sigma2), pnorm(Xn.sorted[i], mean=mu, sd=sigma2) - (i-1)/n))
  }
  return(max(resultat))
}

norme_infinie(data$Age)

#' On peut d'ailleurs la comparer et la confirmer � la fonction ks.test :
ks.test(data$Age, "pnorm", mu, sigma2)$statistic

#' On calcule la valeur de la statique de test
Tn = sqrt(n) * norme_infinie(data$Age)
Tn

#' On simule la loi selon B �chantillons
simulation = matrix(rnorm(B*n, mean=mu, sd=sigma2), n, B)
Ln = sqrt(n) * apply(simulation, 2, norme_infinie)
hist(Ln)
abline(v=Tn, col='red', add=TRUE)

#' Enfin, on calcule la p-value
mean(Tn<Ln)

#' Elle vaut 0.0019, on rejette donc H0. Cette valeur est en accord avec la fonction ks.test
ks.test(data$Age, "pnorm", mu, sigma2)

#' Graphiquement, le "non ajustement" semble v�rifi�
plot(ecdf(data$Age))
curve(pnorm(x, mean=mu, sd=sigma2), add=TRUE)


#' Maintenant qu'on a justifi� le fait de reposer sur des hypoth�ses non-param�triques, on va calculer des intervalles de 
#' confiances et des faire des tests d'hypoth�ses qui utilisent les m�thode bootstrap.  
#'   


#' On note dans ce paragraphe $X_i$ la variable �ge de l'individu $i$. Le param�tre sur lequel on fait l'inf�rence est 
#' $\theta(P) = E_P[X]$. Dans ce cas, on sait que la fonction de perte � condid�rer est : 
#' $$ \gamma_n(\hat{\theta}_n) = \sqrt{n} * (\hat{\theta}_n - \theta) $$ qui est bien une fonction continue de $\hat{\theta}_n$.

#' Commen�ons par le calcul de l'intervalle de confiance pour la moyenne d'�ge des condamn�s.
theta.hat=mean(data$Age)

bootdata.age = boot(data$Age, function(y,i) mean(y[i]), R=B, sim='ordinary')
Ln.hat = sqrt(n) * (bootdata.age$t - theta.hat)

qqplot(Ln.hat, qnorm(seq(0,1,0.01),0,sd(data$Age)))
abline(0,1, col='red')

#' Visuellement on confirme que $\hat{L}_n$ converge bien vers la quantit� pr�vue.

hist(Ln.hat, freq=FALSE)
curve(dnorm(x, sd=sd(data$Age)), add=TRUE, col='blue')

#' L'ajustement est OK. On calcule les quantiles 0.025  et 0.975
quantile(Ln.hat, prob=c(0.025, 0.975)) / sqrt(n)

#' D'o� l'intervalle de confiance final � 95% pour la moyenne d'�ge des condamn�s :
theta.hat + quantile(Ln.hat, prob=c(0.025, 0.975)) / sqrt(n)

#' Passons au subsampling. On commence par �crire la fonction qui cr�e un sous-�chantillon � partir des donn�es :
subsample <- function(data)
{
  data.sampled = sample(data, size=b, replace=FALSE)
}

data.replicated = matrix(data$Age, nrow=n, ncol=B)

data.subsampled = apply(data.replicated, 2, subsample)
theta.subsampled.hat = apply(data.subsampled, 2, mean)

Ln.b.hat = sqrt(b) * (theta.subsampled.hat - theta.hat)

qqplot(Ln.b.hat, qnorm(seq(0,1,0.01), 0, sd(data$Age)))
abline(0,1, col='red')

hist(Ln.b.hat, freq=FALSE)
curve(dnorm(x, sd=sd(data$Age)), add=TRUE, col='blue')

#' On calcule les quantiles 0.025 et 0.975
quantile(Ln.b.hat, prob=c(0.025, 0.975)) / sqrt(b)

#' D'o� l'intervalle de confiance final � 95% pour la moyenne d'�ge des condamn�s :
theta.hat + quantile(Ln.hat, prob=c(0.025, 0.975)) / sqrt(b)

#' Enfin, on compare ce r�sultat � un test de student (dans le cadre param�trique)
t.test(data$Age)

#' Comme pr�vu l'intervalle de confiance est plus petit, mais repose sur une hypoth�se qui n'est dans notre cas pas v�rifi�.  
#' 
#' | Bootstrap | Subsampling | Fisher |
#' |--------------|-------------|-----------|
#' |   38.60     |  35.97 |  38.60       |
#' |--------------|-------------|-----------|
#' |   40.04     |  42.70 |  40.03       |

#' C'est d'ailleurs une remarque g�n�rale que l'on peut faire sur les intervalles de confiance non param�trique : ils sont toujours
#' plus larges que ceux qui sont param�triques, mais poss�dent l'avantage de ne reposer sur aucune hypoth�se de distribution des 
#' donn�es. C'est pourquoi on pr�f�rera faire du non param�trique lorsque la distribution des donn�es est inconnue, et d�s lors
#' que l'on est suffisament sur que les donn�es suivent une certaine distribution, on utilisera les m�thodes param�triques 
#' correspondantes.  
#'   
#'   
#' Enfin, on va tester si la moyenne d'age est de 38 ans en appliquant les deux m�thodes propos�es plus haut.  
#' 
#' Commen�ons par le test bootstrap sous contrainte. On teste : "H0 : $\theta = 38$" contre "H1 : $\theta \neq 38$".
#' La statistique de test est $T_n = |\bar{X_n}|$, et on se fixe $\alpha=0.05$. Afin de se ramener dans le cadre du test "H0 :
#' $\theta = 0$", on va faire "data = data - 38"
#' 

data.38 = data$Age - 38
Tn.hat = sqrt(n) * abs(mean(data.38))

#' On g�n�re l'�chantillon bootstrap sous contrainte, puis on calcule les statistiques de test
data.Age.corrected = data.38 - mean(data.38)
bootdata.age = boot(data.Age.corrected, function(y,i) mean(y[i]), R=B, sim='ordinary')
Tn.star = sqrt(n) * abs(bootdata.age$t)

#' On calcule le quantile empirique qui correspond au niveau $\alpha$ fix�
c_alpha = quantile(Tn.star, prob=1-alpha)

#' On calcule la p-valeur
mean(Tn.star > Tn.hat)

#' Elle est tr�s faible, plus petite que $\alpha$, on rejete donc H0
hist(Tn.star)
abline(v=Tn.hat, col='red')

#' Calculons la puissance du test. Pour ce faire, on va calculer $c_{\alpha}$ pour $\theta = 37, 37.1, ..., 39$ (ici on calcule
#' donc pour $\theta = -1, -0.9, ..., 1$ car on a translat� les donn�es de 38) :
pi = c()
for(theta in seq(-1,1,0.1)){
  data.Age.corrected = data.38 - mean(data.38) + theta
  bootdata.age = boot(data.Age.corrected, function(y,i) mean(y[i]), R=B, sim='ordinary')
  pi = c(pi, mean(sqrt(n) * bootdata.age$t > c_alpha - sqrt(n) * theta) + mean(sqrt(n) * bootdata.age$t < -c_alpha - sqrt(n) * theta))
}
plot(seq(-1,1,0.1), pi, main="Puissance du test")
abline(h=alpha, col='red')

#' On confirme bien que la puissance du test vaut $\alpha$ pour $\theta = 38$, puis augmente et converge vers 1 en s'en �loignant.    
#' Testons par permutation si la moyenne est �gale � 38.
#' Pour ce faire, on va faire data - 38, puis tester "H0 : $\theta = 0$" contre "H1 : $\theta \neq 0$"
data.38 = data$Age - 38
epsilon = rbinom(n, size=1, prob=1/2)
data.38.shuffled = data.38 * epsilon

Tn.hat = sqrt(n) * mean(data$Age - 38)

bootdata.age = boot(data.38.shuffled, function(y,i) mean(y[i]), R=B, sim='ordinary')
Tn.star = sqrt(n) * abs(bootdata.age$t)

#' On calcule la p-valeur
mean(Tn.star > Tn.hat)

#' Elle est aussi tr�s faible, on rejette H0 au niveau 5%
hist(Tn.star, freq=FALSE)
abline(v=Tn.hat, col='red')

#' Etude de la variable jour de la semaine
#' ---
#' 
#' Dans cette derni�re section nous allons nous interesser � la distribution de la variable 'jour de la semaine'. Plus pr�cisement,
#' voici le mod�le que nous allons consid�rer :  
#' Soient $X_1, ..., X_n \sim \mathcal{M}_\theta$ iid les valeurs relev�s des individus dont on dispose. Ils sont donc tir�s
#' depuis une loi multinomiale de param�tre $\theta \in R^7$ avec $|\theta|_1=1$. L'objectif va �tre de pouvoir donner un 
#' intervalle de confiance sur la valeur des composantes de $\theta$ puis de tester si elles sont toutes �gales � $\frac{1}{7}$.

#' Commen�ons par une rapide visualisation. A priori la distribution ne semble pas �tre uniforme.
date = as.Date(data$Date, format='%m/%d/%Y')
date1 = as.factor(sapply(date, weekdays))
plot(date1)

#' Testons le

theta.hat = c(mean(date1 == "lundi"), mean(date1 == "mardi"), mean(date1 == "mercredi"), mean(date1 == "jeudi"), mean(date1 == "vendredi"), mean(date1 == "samedi"), mean(date1 == "dimanche"))

compute_theta <- function(y, i){
  return( c(mean(y[i] == "lundi"), mean(y[i] == "mardi"), mean(y[i] == "mercredi"), mean(y[i] == "jeudi"), mean(y[i] == "vendredi"), mean(y[i] == "samedi"), mean(y[i] == "dimanche")) )
}

bootdata.date1 = boot(date1, compute_theta, R=B, sim='ordinary')

compute_Ln.hat <- function(boot){
  return(sqrt(n) * (boot - theta.hat))
}

Ln.hat = apply(bootdata.date1$t, 1, compute_Ln.hat)

hist(Ln.hat[1,], freq = FALSE)
curve(dnorm(x,sd=sqrt(theta.hat[1]*(1-theta.hat[1]))), add=TRUE, col='blue')

#' Le th�or�me central limite multidimensionnel nous assure que $\sqrt{n} (\hat{\theta}_n - \theta)$ converge vers une loi normale
#' centr�e de matrice de covariance $Diag(\theta_1 (1-\theta_1), \theta_2 (1-\theta_2), ..., \theta_7 (1-`\theta_7))$.

qqplot(Ln.hat[1,], qnorm(seq(0,1,0.01), 0, sqrt(theta.hat[1]*(1-theta.hat[1]))))
abline(0,1, col='red')

#' On peut donc donner un intervalle de confiance sur les param�tres :
theta.hat[1] + quantile(Ln.hat[1,], prob=c(0.025,0.975))/sqrt(n)
theta.hat[2] + quantile(Ln.hat[2,], prob=c(0.025,0.975))/sqrt(n)
theta.hat[3] + quantile(Ln.hat[3,], prob=c(0.025,0.975))/sqrt(n)
theta.hat[4] + quantile(Ln.hat[4,], prob=c(0.025,0.975))/sqrt(n)
theta.hat[5] + quantile(Ln.hat[5,], prob=c(0.025,0.975))/sqrt(n)
theta.hat[6] + quantile(Ln.hat[6,], prob=c(0.025,0.975))/sqrt(n)
theta.hat[7] + quantile(Ln.hat[7,], prob=c(0.025,0.975))/sqrt(n)


#' Terminons avec le test : On veut tester $"H_0 : \theta = (\frac{1}{7}, ..., \frac{1}{7})^T"$ contre 
#' $"H_1 : \theta \neq (\frac{1}{7}, ..., \frac{1}{7})^T"$. On va utiliser le bootstrap sous contrainte.  

Tn.hat = sqrt(n) * abs(theta.hat[1]-1/7)

#' On g�n�re l'�chantillon bootstrap sous contrainte. Ici, vu qu'on suppose connue la classe dans laquelle a �t� g�n�r� les 
#' donn�es (une multinomiale), il suffit de g�n�rer un jeu $X_1^*, ..., X_n^* \sim b(1/7)$ (bernouilli) iid, puis de calculer
#' la valeur de $c_\alpha$ tel que $P[Tn]  puis on calcule les statistiques de test

bootdata.date = rbinom(B, size=n, prob=1/7)
theta.star = bootdata.date/n
Tn.star = sqrt(n) * abs(theta.star-1/7)

#' On calcule le quantile empirique qui correspond au niveau $\alpha$ fix�
c_alpha = quantile(Tn.star, prob=1-alpha)

#' On calcule la p-valeur
mean(Tn.star > Tn.hat)

#' Elle est tr�s faible, plus petite que $\alpha$, on rejete donc H0
hist(Tn.star)
abline(v=Tn.hat, col='red')

#' On pourrait faire ceci pour les 7 jours, tous les tests rejettent $H_0$.  
#' Remarque : on n'a pas appliqu� la proc�dure de Bonferroni ici car le nombre de test � effectuer est n�gligeable devant le
#' nombre d'observations dont on dispose.

