from os import system
import cv2
from TCM_Trabalho_final_classe import *

T = Temperatura2d()
# lado, divisoes
T.Tamanho(1, 30)
# condicao inicial
T.Cond_inicial(0)
# condicao de controno, esquerda, direita, cima, baico
T.Cond_contorno(0, 0, 0, 0)
# funcoes que descrevem o conrtorno,
# inicial caso chame a função Dar_N_passos
# no tempo caso chame a funcao Evolucao
T.Funcoes_contorno(functE=seno, functD=seno)
# partes opcionais como frequencia, inclinacao da reta, etc
T.Opcionais_consotorno(opcE=20, opcD=-20)
T.Difusividade(1)
# limites de plotagem para o eixo z
T.Limite_superior(1)
T.Limite_inferior(-1)
# cores disponiveis:
# 0: azul -> vermelho
# 1: verde -> amarelo
# 2: rubro -> branco
T.Cor(2)

# frames por segundo
fps = 10.0
# passos por frame
passo = 15
total_frames = 1000
# não alterar !!
unidade_de_passo = T.dt
dimensao = (640, 480)
nome_video = 'Video_TCM_tab_final.avi'

video = cv2.VideoWriter(nome_video, cv2.VideoWriter_fourcc(*'mp4v'), fps, dimensao)

for t in range(0, total_frames):
    tempo = t*unidade_de_passo
    T.Salvar('grafico{}.jpg'.format(''))
    frame = cv2.imread('grafico{}.jpg'.format(''))
    frame = cv2.resize(frame, dsize=dimensao, interpolation=cv2.INTER_CUBIC)
    video.write(frame)
    # funcao a ser feita por loop
    T.Evolucao(passo, fe=1, fc=1, fb=1, fd=1)
    system("cls")
    print("Fazendo video...")
    print(str(100*t/total_frames) + "%")
print("Video concluido!")

video.release()
cv2.destroyAllWindows()
