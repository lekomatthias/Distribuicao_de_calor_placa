import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from numba import njit

def const(x, opc=0):
    return 0

def constx(x, opc=0):
    return 1*opc

def reta(t, opc=1):
    # opc == inclinacao
    return opc*t

def seno(pos, opc=1):
    # opc == frequencia
    return np.sin(pos*opc)

# coloca a condicao inicial da placa
def cond_inic(vetor2d, valor, dx, functX=const, functY=const):
    tamanho = len(vetor2d)
    for j in range(0, tamanho):
        for i in range(0, tamanho):
             vetor2d[j, i] = valor + functX(i*dx) + functY(j*dx)
    return vetor2d

# coloca as condicoes de contorno na placa pelo tempo
def cond_cont(vetor2d, valorE, valorD, valorC, valorB, lado=1, 
                functE=const, functD=const, functC=const, functB=const,
                opcE=1, opcD=1, opcC=1, opcB=1):
    final = len(vetor2d)
    variacao = lado/final
    for i in range(0, final):
        vetor2d[i, 0] = valorB + functB(i*variacao, opcB)
    for i in range(1, final):
        vetor2d[i, final-1] = valorC + functC(i*variacao, opcC)
    for i in range(1, final):
        vetor2d[0, i] = valorE + functE(i*variacao, opcE)
    for i in range(1, final):
        vetor2d[final-1, i] = valorD + functD(i*variacao, opcD)
    return vetor2d

# usa o método de diferenças finitas fixandos as bordas
# para decidir a temperatura em cada ponto do tempo seguinte
@njit
def evolucao_passo(Temperatura, dx, dt, difus):
    aux = Temperatura.copy()
    N = len(Temperatura)
    for j in range (1, N - 1):
        for i in range (1, N - 1) :
            # influencia horizontal
            h = (dt/(dx*dx))*(Temperatura[i+1, j] - 2*Temperatura[i, j] + Temperatura[i-1, j])
            # influencia vertical
            v = (dt/(dx*dx))*(Temperatura[i, j+1] - 2*Temperatura[i, j] + Temperatura[i, j-1])
            # derivada segunda pela equacao do calor 2d
            aux[i, j] = Temperatura[i, j] + difus*(h + v)
    return aux

# faz a evoluçao complete para um tempo definido
@njit
def evolucao_geral(Temperatura, dx, dt, dif, t):
    for _ in range(0, t+1):
        Temperatura = evolucao_passo(Temperatura, dx, dt, dif)
    return Temperatura

def dirichlet_esquerda(Temperatura, valor, dt, t_atual, final, funct=const, opc=1):
    for pos in range(0, final):
        Temperatura[0, pos] = valor + funct(dt*t_atual*opc)
    return Temperatura

def dirichlet_direita(Temperatura, valor, dt, t_atual, final, funct=const, opc=1):
    for pos in range(0, final):
        Temperatura[final, pos] = valor + funct(dt*t_atual*opc)
    return Temperatura

def dirichlet_cima(Temperatura, valor, dt, t_atual, final, funct=const, opc=1):
    for pos in range(0, final):
        Temperatura[pos, final] = valor + funct(dt*t_atual*opc)
    return Temperatura

def dirichlet_baixo(Temperatura, valor, dt, t_atual, final, funct=const, opc=1):
    for pos in range(0, final):
        Temperatura[pos, 0] = valor + funct(dt*t_atual*opc)
    return Temperatura

def cond_cont_variada(vetor2d, valorE, valorD, valorC, valorB, dt, tempo, 
                    functE=const, functD=const, functC=const, functB=const,
                    opcE=1, opcD=1, opcC=1, opcB=1):
    final = len(vetor2d)-1
    dirichlet_esquerda(vetor2d, valorE, dt, tempo, final, functE, opcE)
    dirichlet_direita(vetor2d, valorD, dt, tempo, final, functD, opcD)
    dirichlet_cima(vetor2d, valorC, dt, tempo, final, functC, opcC)
    dirichlet_baixo(vetor2d, valorB, dt, tempo, final, functB, opcB)
    return vetor2d

def evolucao_dirichlet_variada(Temperatura, dx, dt, dif, t_atual, t_final,
                                tempiE=0, tempiD=0, tempiC=0, tempiB=0,
                                varE=const, varD=const, varC=const, varB=const,
                                opcE=1, opcD=1, opcC=1, opcB=1):
    for i in range(t_atual, t_final+1):
        cond_cont_variada(Temperatura, tempiE, tempiD, tempiC, tempiB, dt, i*dt,
                            functB=varB, functC=varC, functD=varD, functE=varE,
                            opcB=opcB, opcC=opcC, opcD=opcD, opcE=opcE)
        Temperatura = evolucao_passo(Temperatura, dx, dt, dif)
    return Temperatura

def neumann_esquerda(Temperatura, dx, dt, t_atual, final, fluxo=const, opc=1):
    for pos in range (0, final):
        Temperatura[0, pos] = fluxo(t_atual*dt, opc=opc)
        Temperatura[0, pos] = -2*(Temperatura[0, pos]*dx + Temperatura[2, pos]/2 - 2*Temperatura[1, pos])/3
    return Temperatura

def neumann_direita(Temperatura, dx, dt, t_atual, final, fluxo=const, opc=1):
    for pos in range (0, final):
        Temperatura[final, pos] = fluxo(t_atual*dt, opc=opc)
        Temperatura[final, pos] = -2*(Temperatura[final, pos]*dx + Temperatura[final-2, pos]/2 - 2*Temperatura[final-1, pos])/3
    return Temperatura
        
def neumann_cima(Temperatura, dx, dt, t_atual, final, fluxo=const, opc=1):
    for pos in range (0, final):
        Temperatura[pos, final] = fluxo(t_atual*dt, opc=opc)
        Temperatura[pos, final] = -2*(Temperatura[pos, final]*dx + Temperatura[pos, final-2]/2 - 2*Temperatura[pos, final-1])/3
    return Temperatura

def neumann_baixo(Temperatura, dx, dt, t_atual, final, fluxo=const, opc=1):
    for pos in range (0, final):
        Temperatura[pos, 0] = fluxo(t_atual*dt, opc=opc)
        Temperatura[pos, 0] = -2*(Temperatura[pos, 0]*dx + Temperatura[pos, 2]/2 - 2*Temperatura[pos, 1])/3
    return Temperatura

def evolucao_neumann(Temperatura, dx, dt, dif, t_atual, t_final,
                    fluxoE=const, fluxoD=const, fluxoC=const, fluxoB=const,
                    opcE=1, opcD=1, opcC=1, opcB=1):
    final = len(Temperatura)-1
    for k in range (t_atual, t_final+1):
        # derivadas nas bordas para que ela o fluxo seja o dado na borda
        neumann_baixo(Temperatura, dx, dt, k, final, fluxo=fluxoB, opc=opcB)
        neumann_cima(Temperatura, dx, dt, k, final, fluxo=fluxoC, opc=opcC)
        neumann_esquerda(Temperatura, dx, dt, k, final, fluxo=fluxoE, opc=opcE)
        neumann_direita(Temperatura, dx, dt, k, final, fluxo=fluxoD, opc=opcD)
        
        Temperatura = evolucao_passo(Temperatura, dx, dt, dif)

    return Temperatura

# plota o gráfico da temperatura e usa os outros dados para
# colocar legendas adequadas nele
def plotTemp3d(Temperatura, L, t, limite=[], cor=1):
    # criacao do espaco 3d
    X, Y = np.mgrid[0:L:(len(Temperatura[0]))*1j, 0:L:(len(Temperatura))*1j]
    try:
        plt.style.use ('default')
    except:
        pass
    ax = plt.axes(projection = '3d')

    if cor == 0:
        ax.plot_surface(X, Y, Temperatura, cmap = cm.coolwarm)
    elif cor == 1:
        ax.plot_surface(X, Y, Temperatura, cmap = 'viridis')
    else:
        ax.plot_surface(X, Y, Temperatura, cmap = 'hot')

    ax.set_title(' lado: %.3fm '%L
                +'\ntempo: %.3fs '%t)
    ax.set_ylabel('$y$')
    ax.set_xlabel('$x$')
    ax.set_zlabel('$T[°C]$')
    if limite[0] != limite[1]:
        ax.set_zlim(limite[0],limite[1])

    plt.show()
    # plt.cla()

# @jit(nopython=True)
def Salvar_figura(Temperatura, Largura, tempo, dt, nome, limite=[], cor=1):
    # criacao do espaco 3d
    X, Y = np.mgrid[0:Largura:(len(Temperatura[0]))*1j, 0:Largura:(len(Temperatura))*1j]

    fig = plt.figure() # tentativa de nao criar figura e nao correr o risco de nao deletar
    ax = plt.axes(projection = '3d')
    if cor == 0:
        ax.plot_surface(X, Y, Temperatura, cmap = cm.coolwarm)
    elif cor == 1:
        ax.plot_surface(X, Y, Temperatura, cmap = 'viridis')
    else:
        ax.plot_surface(X, Y, Temperatura, cmap = 'hot')
    ax.set_title(' lado: %.3fm '%Largura
                +'\ntempo: %.3fs '%(tempo*dt))
    ax.set_ylabel('$y$')
    ax.set_xlabel('$x$')
    ax.set_zlabel('$T[°C]$')
    if limite[0] != limite[1]:
        ax.set_zlim(limite[0],limite[1])

    plt.savefig(nome, dpi=200, bbox_inches='tight')
    plt.cla()
    plt.close(fig) # garantia de que deletei  a figura
