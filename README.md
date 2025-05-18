<div align="center"><a href="https://github.com/IsaacAlves7/cpp-programming"><img src="https://user-images.githubusercontent.com/61624336/181665378-52c87f52-fdc1-4ca9-b4bf-e8f9cd523d5d.png"></a></div>

# It's a repository of C/C++ languages 🅲

> 🅲 **Preparação**: Para este conteúdo, o aluno deverá dispor de um computador com acesso à internet, um web browser com suporte a HTML 5 (Google Chrome, Mozilla Firefox, Microsoft Edge, Safari, Opera etc.), um editor de texto ou IDE (VSCode etc.) e o software C/C++, com a versão mais recente, instalado na sua máquina local.

- https://exercism.org/tracks/c
- http://www.dpi.inpe.br/~carlos/Academicos/Cursos/LinguagemC/Cap_1.html
- https://exercism.org/tracks/cpp
- https://cplusplus.com/

# 🐒 Linguagem de programação

<div align="center"><img src="https://user-images.githubusercontent.com/61624336/112900537-065ce480-90ba-11eb-86f7-f9006445876a.png"></div>

Hoje em dia, o desenvolvimento de sistemas se baseia em vários e diferentes paradigmas, tais como os listados a seguir:

- **Imperativo (Procedural)**: Segue sequências de comandos ordenados segundo uma lógica.
- **Funcional**: Trabalha com a divisão de problemas através de funções, que resolvem separadamente problemas menores e que, ao serem organizados, resolvem o problema como um todo.
- **Lógico**: Voltado ao desenvolvimento de problemas de lógica e usado em sistemas de inteligência computacional.
- **Orientado a Objetos (OO)**: Define um conjunto de classes para dividir o problema e realiza a interação entre as diferentes classes para também resolver o problema como um todo.

# 🅲 linguagem C
<img src="https://upload.wikimedia.org/wikipedia/commons/1/18/C_Programming_Language.svg" height="77" align="right">

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

Exemplo:
Podemos detalhar o processo exemplificando a compilação em Linux de um programa em C simples como o abaixo, de nome `compilacao.c`, que escreve na tela a frase `Hello, World!`:
```c
#include <stdio.h>
int main()
{
printf("Hello World!");
    return 0;
}
```

O comando `#include` serve para incluir uma biblioteca e o comando `<stdio.h>` serve para a entrada e saída de dados;
A função `main ()` e `{}` dentro dela o comando `printf ("Ola mundo \n");`  
Necessita do `;` para rodar a função!

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

### Pré-processador
O **pré-processador** é um pequeno software que aceita o arquivo-fonte C e executa as tarefas abaixo.

- Remove comentários do código-fonte.
- Faz a expansão dos arquivos de cabeçalho incluídos.
- Gera um arquivo temporário com a extensão `.i` após o pré-processamento. Ele insere o conteúdo dos arquivos de cabeçalho no arquivo de código-fonte. O arquivo gerado pelo pré-processador é maior do que o arquivo de origem original.

### Compilador
Na próxima fase da compilação C, o compilador entra em ação. Ele aceita o arquivo pré-processado temporário nome_do_arquivo.i gerado pelo pré-processador e executa as seguintes tarefas:

- Verifica o programa C para erros de sintaxe.
- Traduz o arquivo em código intermediário, ou seja, em linguagem assembly.
- Otimiza,opcionalmente, o código traduzido para melhor desempenho.
- Gera um código intermediário na linguagem assembly, após a compilação, como `nome_do_arquivo.s`. É a versão de montagem do código-fonte.

### Montador (Assembler)
Passando para a próxima fase de compilação, o **assembler** aceita o código-fonte compilado (`nome_do_arquivo.s`) e o traduz em código de máquina de baixo nível. Após a montagem bem-sucedida, gera o arquivo `nome_do_arquivo.o` (no Linux) ou `nome_do_arquivo.obj` (no Windows) conhecido como **arquivo objeto**. No nosso caso, gera o arquivo **compilacao.o**.

### Vinculador (Linker)
Finalmente, o **linker** entra em ação e executa a tarefa final do processo de compilação. Aceita o arquivo intermediário `nome_do_arquivo.o` gerado pelo assembler.

Ele liga todas as chamadas de função com sua definição original. O que significa que a função `printf ()` é vinculada à sua definição original. O **vinculador** gera o arquivo executável final.

## Variáveis e tipos de dados
Na programação, **uma variável** é **um contêiner** (área de armazenamento) para armazenar dados.

Para indicar a área de armazenamento, cada variável deve receber um nome exclusivo (identificador). Os **nomes de variáveis** são apenas a representação simbólica de um local de memória.

Exemplo:

```c
int resultado = 95;
```

Aqui, resultado é uma variável do tipo **inteiro** (`int`). Para esta variável, é atribuído um valor inteiro, `95`.

O valor de uma variável pode ser alterado, como abaixo. Daí o nome, **variável**.

```c
char ch = 'a';

// algum código

ch = 'l';
```

## Regras para nomear uma variável
Um nome de variável pode ter letras ( A primeira letra de uma variável deve ser uma letra), dígitos e símbolo "_". 

<blockquote><b>ATENÇÃO!</b> Não há nenhuma regra sobre o tamanho que um nome de variável (identificador) pode ter. No entanto, podemos ter problemas em alguns compiladores se o nome da variável tiver mais de 31 caracteres.</blockquote>

**C** é uma linguagem fortemente tipada ou tipificada. Isso significa que o tipo da variável não pode ser alterado depois de declarado.

Exemplo:
```c
intnumero  = 5;                                // variável inteira

numero = 5.5;                                   // erro

floatnumero ;                                    // erro
```

Aqui, o tipo de variável numérica é `int`. Você não pode atribuir um valor de **ponto flutuante** (`5.5`) a essa variável. Além disso, você não pode redefinir o tipo da variável para `float`.

<blockquote>A propósito, para armazenar valores com casas decimais em C, você precisa declarar seu tipo para <code>double</code> ou <code>float</code>.</blockquote>

## Constantes 
Uma **constante** é um valor (ou um identificador) cujo valor não pode ser alterado em um programa.

<blockquote>
1, 2.5, 'c' etc.

Aqui, 1, 2.5 e 'c' são constantes literais. Não se pode atribuir valores diferentes a esses termos.

constfloat PI = 3,14;

Observe que adicionamos a palavra-chave const.

Aqui, PI é uma constante simbólica. Na verdade, é uma variável, no entanto, seu valor não pode ser alterado.</blockquote>

Tipos de constantes:

Veja os tipos de constantes que podem ser usadas em C:

- Constantes inteiras.
- Constantes de ponto flutuante.
- Constantes de caracteres.

Uma constante de caractere é criada, colocando-se um único caractere entre aspas simples.

Por exemplo: 'a', 'm', 'F', '2', '}' etc.

- Sequências de escape

Às vezes, é necessário usar caracteres que não podem ser digitados ou que tenham significado especial na programação C. Para usar esses caracteres, a sequência de escape é usada.

Por exemplo: `\n` é usado para nova linha. `\t` como tabulação horizontal. A barra invertida (`\`) faz com que se escape do modo normal, em que os caracteres são manipulados pelo compilador.

- String literal
Uma **string literal** é uma sequência de caracteres entre aspas duplas.

Exemplo:

```c

"legal"          // constante de string

""               // constante de cadeia nula

"      "             // constante de seis espaços em branco

"A"              // constante de string com caractere único

"Resultado eh\n"     // imprime string com nova linha
```

- Enumerações
A palavra-chave `enum` é usada para definir tipos de enumeração.

Exemplos:
```c
enum cor {amarelo, verde, preto, branco};
```

<blockquote>Aqui, a cor é uma variável e amarelo, verde, preto e branco são as constantes de enumeração com valor 0, 1, 2 e 3, respectivamente.</blockquote>

Pode-se definir constantes simbólicas usando-se também a palavra **#define**.

## Tipos de dados e modificadores
São 5 os tipos de dados básicos em C:

<table>
    <tr>
      <td><code>char</code></td>
      <td>Caractere</td>
    </tr>
    <tr>
        <td><code></code></td>
    </tr>
</table>

# [C] `Hello, World`
```c
# include <stdio.h>

int main(void){
  printf("Hello, World!\n");
}
```

## Compilação do código em C
**Entrada (Input)**:
```c
make main.c
./main
```

**Saída (Input)**:
```c

Hello, World!
```

# Linguagem C++
<img src="https://cdn.worldvectorlogo.com/logos/c.svg" height="77" align="right">

A linguagem C++ (em português: Pronuncia-se "cê mais mais") é uma linguagem de programação compilada multi-paradigma (seu suporte inclui linguagem imperativa, orientada a objetos e genérica) e de uso geral. Desde os anos 1990 é uma das linguagens comerciais mais populares, sendo bastante usada também na academia por seu grande desempenho e base de utilizadores.

Bjarne Stroustrup desenvolveu o C++ (originalmente com o nome C with Classes, que significa C com classes em português) em 1983 no Bell Labs como um adicional à linguagem C. Novas características foram adicionadas com o tempo, como funções virtuais, sobrecarga de operadores, herança múltipla, gabaritos e tratamento de exceções. Após a padronização ISO realizada em 1998 e a posterior revisão realizada em 2003, uma nova versão da especificação da linguagem foi lançada em dezembro de 2014, conhecida informalmente como C++14.

C++ é uma linguagem de programação de propósito geral, amplamente utilizada em diversas áreas como desenvolvimento de software, jogos, sistemas embarcados e inteligência artificial. É uma linguagem compilada, de alto nível e que permite um controle mais direto sobre o hardware, tornando-a muito eficiente. 

## [C++] `Hello, World!`


# [C++] Sistemas Digitais
As áreas de TI e Comunicação trazem, a todo o momento, modificações, inovações, adequações, enfim, apresentam-se de forma cada vez mais interessantes para o usuário e
desafiadoras para o profissional que as constrói. Assim, preparar equipes capazes de conceber, planejar e desenvolver soluções que funcionarão nas futuras gerações das áreas
de Tecnologia da Informação e Comunicação (TIC) apresenta-se como demanda urgente aos Cursos da área de TI e um desafio às práticas pedagógicas do professor para o ensino
da computação.

Circuitos digitais são definidos como circuitos eletrônicos que empregam a utilização de sinais elétricos em apenas dois níveis de corrente (ou tensão) para definir a representação de valores binários. A importância do estudo dos circuitos lógicos como base para o estudo dos sistemas digitais é de grande relevância, uma vez que são a base dos circuitos encontrados nos computadores atuais e em uma enorme quantidade de dispositivos e instrumentos usados em todas as áreas.

# [C++] Programação de Microcontroladores
<img src="https://upload.wikimedia.org/wikipedia/commons/5/58/DIL16_Labelled.svg" height="77" align="right">

A **programação de microcontroladores** envolve escrever códigos que permitem que esses pequenos computadores internos controlem dispositivos eletrônicos, desde sistemas industriais até dispositivos de uso diário como smartphones. Essa programação é fundamental para o funcionamento de diversas aplicações e sistemas, sendo essencial para a criação de soluções em áreas como automação, eletrônica e IoT. 

# Arduíno
<img src="https://upload.wikimedia.org/wikipedia/commons/8/87/Arduino_Logo.svg" height="77" align="right">

# QT
<img src="https://cdn.worldvectorlogo.com/logos/qt-1.svg" height="77" align="right">

O **QT** é um conjunto de ferramentas multiplataforma que combina uma estrutura de desenvolvimento de software e uma interface gráfica (GUI) altamente eficiente, amplamente utilizado para criar aplicativos com interfaces de usuário modernas e responsivas. Originalmente desenvolvido pela empresa norueguesa Trolltech, agora mantido pela empresa The Qt Company, o QT se destaca por sua flexibilidade, desempenho e pela facilidade de desenvolvimento de aplicações que podem ser executadas em diversos sistemas operacionais, como Windows, macOS, Linux, Android e iOS, sem necessidade de reescrita do código. Em resumo, o QT é uma solução robusta e versátil para o desenvolvimento de software, capacitando desenvolvedores a criar aplicativos potentes, esteticamente agradáveis e executáveis em várias plataformas com um único código-base. Ele continua a ser uma escolha popular em indústrias que exigem alta performance, flexibilidade e consistência entre dispositivos e sistemas operacionais.

Um dos principais atrativos do QT é sua abordagem orientada a objetos, que é implementada principalmente em C++. Ele fornece um conjunto abrangente de bibliotecas e APIs que simplificam o desenvolvimento de componentes gráficos, manipulação de eventos, acesso a banco de dados, redes e até mesmo processamento multithreading. Além disso, a integração de seu mecanismo de sinais e slots facilita a comunicação entre diferentes componentes do software, tornando o código mais modular e organizado.

Outro aspecto notável do QT é sua compatibilidade com linguagens de script, como Python, por meio do PyQt ou PySide, o que o torna acessível a desenvolvedores que preferem uma curva de aprendizado mais suave. Sua documentação rica e comunidade ativa também são grandes vantagens para iniciantes e veteranos que precisam de suporte ou recursos adicionais para seus projetos.

O QT também brilha no desenvolvimento de aplicações embarcadas, especialmente em dispositivos IoT, automação industrial e sistemas automotivos, devido ao seu desempenho leve e sua capacidade de criar interfaces gráficas de alta qualidade mesmo em hardware limitado. Suas ferramentas incluem o Qt Creator, um ambiente de desenvolvimento integrado (IDE) que simplifica o design visual de interfaces, a prototipagem e a depuração.

Embora seja amplamente elogiado, o QT também apresenta desafios. O licenciamento comercial pode ser caro, especialmente para pequenas empresas ou desenvolvedores individuais, embora exista uma versão de código aberto disponível sob a licença LGPL. A complexidade inicial para quem não está acostumado com C++ ou com frameworks avançados também pode ser um obstáculo.
