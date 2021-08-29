# <img src="https://cdn.worldvectorlogo.com/logos/c.svg" height="27"> It's a repository of C/ C++ programming 🅲🔢

<blockquote>This repository contains Full-Stack development in C/ C++ languages!</blockquote>

<div align="center"><img src="https://res.cloudinary.com/practicaldev/image/fetch/s--3u1aWUCM--/c_imagga_scale,f_auto,fl_progressive,h_420,q_auto,w_1000/https://dev-to-uploads.s3.amazonaws.com/i/stfvlecgmmp4dso3v0iv.jpg"></div>

# 🅲 A História da linguagem C 🅲
<div align="center"><img src="https://upload.wikimedia.org/wikipedia/commons/1/18/C_Programming_Language.svg" height="177"></div><br \>

Uma tarefa difícil na área de computação é convencer um estudante que aprender uma nova linguagem de programação, ou usar uma linguagem que não é a preferida dele, é necessário e essencial dentro de uma disciplina. Quando se trata de uma linguagem que para alguns está ultrapassada, como a **linguagem C**, a tarefa é ainda mais difícil.

Existem muitas razões para o aprendizado de C ser fundamental. Muitos a consideram até a "mãe de todas as linguagens de programação".

Ela foi projetada para implementar o Sistema Operacional **Unix**, ficando próxima ao sistema operacional, o que a torna uma linguagem eficiente devido ao seu hábil gerenciamento de recursos no nível do sistema.

Outro ponto importante é que essa linguagem não é limitada, mas amplamente utilizada em:

- Sistemas operacionais.
- Compiladores de linguagem.
- Drivers de rede.
- Interpretadores de linguagem.
- Áreas de desenvolvimento de utilitários de sistema.
- Sistemas embarcados (embutidos).

Outras vantagens da linguagem C, incluem:

- **Onipresente**: Qualquer que seja a plataforma, C provavelmente está disponível.
- **Portável**: Um programa em C compila com modificações mínimas em outras plataformas − às vezes até funciona de imediato.
- **Simples**: C é muito simples de aprender e praticamente não requer dependências. Basta um simples PC com o compilador e tudo está pronto para criar programas.

### Estrutura de um programa em C e processo de compilação
**C** é uma linguagem considerada de nível intermediário e precisa de um compilador para criar um código executável e para que o programa possa funcionar em uma máquina.

<blockquote><b>Compilação</b> é processo de tradução do código-fonte escrito para um código de máquina. É feita por um software especial conhecido como <b>compilador</b>, que verifica o código-fonte em busca de qualquer erro sintático ou estrutural e gera um código-objeto com extensão <code>.obj</code> (no Windows) ou <code>.o</code> (no Linux), se o código-fonte estiver livre de erros.

Iremos utilizar um compilador para o Windows, o <a href="">MinGW</a>.</blockquote>

Toda a compilação é dividida em quatro etapas:

1. **Pré-processamento**.
2. **Compilação**.
3. **Montagem (assembler)**.
4. **Vinculação (linker)**.

A figura 1 descreve todo o processo de compilação em C.

<div align="center"><img src="https://estacio.webaula.com.br/cursos/go0374/galeria/aula1/img/figura1.svg"></div>

<h5 align="center"><i>Figura 1: Processo de compilação de um programa</i></h5>

Uma **IDE**, ou **Ambiente de Desenvolvimento Integrado** (Integrated Development Environment), reúne características e ferramentas de apoio ao desenvolvimento de software com o objetivo de agilizar este processo, disponibilizando todo o processo de compilação no apertar de um botão.

### Exemplo:
Podemos detalhar o processo exemplificando a compilação em Linux de um programa em C simples como o abaixo, de nome `compilacao.c`, que escreve na tela a frase `Hello, World!`:
```c
#include <stdio.h>
int main()
{
printf("Hello World!");
    return 0;}
```
Para compilar o programa acima abre-se o prompt de comando e pressiona-se o comando abaixo:
```
gcc -save-tempscompilacao.c -o compilacao
```

A opção `-save-temps` preservará e salvará todos os arquivos temporários criados durante a compilação em **C**. Ele gerará quatro arquivos no mesmo diretório:

- compilacao.i (gerado pelo pré-processador).
- compilacao.s (gerado pelo compilador).
- compilacao.o (gerado pelo montador).
- compilacao (no Linux gerado pelo linker) ou (compilacao.exe no Windows)

Agora, entenda o papel de cada elemento do processo de compilaÇão:


## linguagem C++ (CPP - CPlusPlus)
<div align="center"><img src="https://cdn.worldvectorlogo.com/logos/c.svg" height="177"></div><br \>

## linguagem C# (C-Sharp)
<div align="center"><img src="https://iconape.com/wp-content/files/sh/51404/svg/c--4.svg" height="177"></div><br \>
