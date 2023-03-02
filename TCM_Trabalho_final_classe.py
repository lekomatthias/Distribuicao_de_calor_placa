from TCM_Trabalho_final_lib import *

class Temperatura2d():

    def __init__(self):

        self.tempo = 0
        self.L = 1
        self.C = self.L
        self.N = 3
        self.Ti = 0
        # esquerda, direita, cima, baixo
        self.contornos = [0, 0, 0, 0]
        self.funct_contorno = [const, const, const, const]
        self.opcionais_contorno = [1, 1, 1, 1]
        self.difusividade = 0
        self.cor = 0
        self.dx = self.L/self.N
        self.dt = 0.1*self.dx*self.dx
        self.limites = [0, 0]
        self.temperatura = np.zeros((self.N+1, self.N+1), dtype=float)

    def Tamanho(self, tamanho, divisoes):
        self.L = tamanho
        self.C = tamanho
        self.N = divisoes+1
        self.dx = self.L/self.N
        self.dt = 0.1*self.dx*self.dx
        self.temperatura = np.zeros((self.N+1, self.N+1), dtype=float)

    def Cond_inicial(self, temperatura, functX=const, functY=const):
        self.Ti = temperatura
        self.temperatura = cond_inic(self.temperatura, temperatura, self.dx, functX=functX, functY=functY)

    def Cond_contorno(self, esquerda, direita, cima, baixo, lado=1):
        self.contornos = [esquerda, direita, cima, baixo]
        self.temperatura = cond_cont(self.temperatura, esquerda, direita, cima, baixo, lado=lado,
                                    functE=self.funct_contorno[0], functD=self.funct_contorno[1],
                                    functC=self.funct_contorno[2], functB=self.funct_contorno[3],
                                    opcE=self.opcionais_contorno[0], opcD=self.opcionais_contorno[1],
                                    opcC=self.opcionais_contorno[2], opcB=self.opcionais_contorno[3])
    
    def Funcoes_contorno(self, functE=const, functD=const, functC=const, functB=const):
        self.funct_contorno = [functE, functD, functC, functB]

    def Opcionais_consotorno(self, opcE=1, opcD=1, opcC=1, opcB=1):
        self.opcionais_contorno = [opcE, opcD, opcC, opcB]

    def Difusividade(self, difusividade):
        self.difusividade = difusividade

    def Limite_superior(self, limite):
        try:
            limite_num = float(limite)
        except:
            limite_num = self.limites[1]
        # coloca a maior temperatura como limite superior
        if limite == '':
            self.limites[1] = max(self.contornos)
            if self.Ti > self.limites[1]:
                self.limites[1] = self.Ti
        else:
            self.limites[1] = limite_num

    def Limite_inferior(self, limite):
        try:
            limite_num = float(limite)
        except:
            limite_num = self.limites[0]
        # coloca a menor temperatura como limite inferior
        if limite == '':
            self.limites[0] = min(self.contornos)
            if self.Ti < self.limites[0]:
                self.limites[0] = self.Ti
        else:
            self.limites[0] = limite_num

    def Ajuste_cor(self, valorsup, valorinf):
        final = len(self.temperatura)-1
        self.temperatura[0, final] = valorsup
        self.temperatura[final, final] = valorinf

    def Cor(self, cor):
        self.cor = cor

    def Passo(self):
        return self.dt

    def Dar_passo(self):
        self.temperatura = evolucao_passo(self.temperatura, self.dx, self.dt, self.difusividade)
        self.tempo = self.tempo + 1
        return self.temperatura

    def Dar_N_passos(self, numero):
        self.temperatura = evolucao_geral(self.temperatura, self.dx, self.dt, self.difusividade, numero)
        self.tempo = self.tempo + numero
        return self.temperatura

    def Evolucao(self, numero, fe=0, fd=0, fc=0, fb=0):
        for _ in range(0, numero):
            if fe == 1:
                self.temperatura = neumann_esquerda(self.temperatura, self.dx, self.dt,
                                                self.tempo, self.N, self.funct_contorno[0], self.opcionais_contorno[0])
            else:
                self.temperatura = dirichlet_esquerda(self.temperatura, self.contornos[0], self.dt,
                                                    self.tempo, self.N, self.funct_contorno[0], self.opcionais_contorno[0])
            if fd == 1:
                self.temperatura = neumann_direita(self.temperatura, self.dx, self.dt,
                                                self.tempo, self.N, self.funct_contorno[1], self.opcionais_contorno[1])
            else:
                self.temperatura = dirichlet_direita(self.temperatura, self.contornos[1], self.dt,
                                                    self.tempo, self.N, self.funct_contorno[1], self.opcionais_contorno[1])
            if fc == 1:
                self.temperatura = neumann_cima(self.temperatura, self.dx, self.dt,
                                                self.tempo, self.N, self.funct_contorno[2], self.opcionais_contorno[2])
            else:
                self.temperatura = dirichlet_cima(self.temperatura, self.contornos[2], self.dt,
                                                    self.tempo, self.N, self.funct_contorno[2], self.opcionais_contorno[2])
            if fb == 1:
                self.temperatura = neumann_baixo(self.temperatura, self.dx, self.dt,
                                                self.tempo, self.N, self.funct_contorno[3], self.opcionais_contorno[3])
            else:
                self.temperatura = dirichlet_baixo(self.temperatura, self.contornos[3], self.dt,
                                                    self.tempo, self.N, self.funct_contorno[3], self.opcionais_contorno[3])
            self.temperatura = evolucao_passo(self.temperatura, self.dx, self.dt, self.difusividade)
            self.tempo = self.tempo + 1
        return self.temperatura

    def Dar_N_passos_variando_borda(self, numero, varE=const, varD=const, varC=const, varB=const,
                                                    opcE=1, opcD=1, opcC=1, opcB=1):
            self.temperatura = evolucao_dirichlet_variada(self.temperatura, self.dx, self.dt, self.difusividade,
                                                            self.tempo, self.tempo+numero,
                                                            self.contornos[0], self.contornos[1], self.contornos[2], self.contornos[3],
                                                            varB=varB, varC=varC, varD=varD, varE=varE,
                                                            opcB=opcB, opcC=opcC, opcD=opcD, opcE=opcE)
            self.tempo = self.tempo + numero
            return self.temperatura

    def Dar_N_passos_Neumann(self, numero, fluxoE=const, fluxoD=const, fluxoC=const, fluxoB=const,
                                            opcE=1, opcD=1, opcC=1, opcB=1):
        self.temperatura = evolucao_neumann(self.temperatura, self.dx, self.dt, 
                                            self.difusividade, self.tempo, self.tempo + numero, 
                                            fluxoC=fluxoC, fluxoB=fluxoB, fluxoD=fluxoD, fluxoE=fluxoE,
                                            opcE=opcE, opcD=opcD, opcC=opcC, opcB=opcB)
        self.tempo = self.tempo + numero
        return self.temperatura

    def Salvar(self, nome):
        Salvar_figura(self.temperatura, self.L, self.tempo, self.dt,nome, limite=self.limites, cor=self.cor)

    def plotar(self):
        plotTemp3d(self.temperatura, self.L, self.tempo*self.dt, self.limites, self.cor)