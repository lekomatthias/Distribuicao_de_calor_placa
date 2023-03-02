from TCM_Trabalho_final_classe import *

# largura da barra
L = 1
C = L
# divisoes
N = 30
# amostras no tempo
t = 750
# condicoes iniciais
Ti = 0
# condicoes de contorno esquerda, direita, cima, baixo
# e suas respectivas funcoes
Te = 0
Fe = seno
opce = 20
Td = 0
Fd = seno
opcd = -20
Tc = 0
Fc = const
opcc = 1
Tb = 0
Fb = const
opcb = 1
# difusividade do material
difu = 1
cor = 2

T = Temperatura2d()
T.Tamanho(L, N)
T.Cond_inicial(Ti)
T.Cond_contorno(Te, Td, Tc, Tb, lado=L)
T.Funcoes_contorno(functE=Fe, functD=Fd, functC=Fc, functB=Fb)
T.Opcionais_consotorno(opcE=opce, opcD=opcd, opcC=opcc, opcB=opcb)
T.Difusividade(difu)
T.Cor(cor)
T.Evolucao(t,fe=1 , fd=1, fc=1, fb=1)
# T.Dar_N_passos(t)
T.Limite_superior(1)
T.Limite_inferior(-1)
# T.Ajuste_cor(2, 0)
T.plotar()
# T.Salvar('pp4-2.pdf')
