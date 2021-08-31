# <img src="https://upload.wikimedia.org/wikipedia/commons/1/18/C_Programming_Language.svg" height="27"> It's a repository of C programming 🅲🔢

<blockquote>This repository contains Full-Stack development in C/ C++ languages!</blockquote>

<div align="center"><img src="https://askatul.com/wp-content/uploads/2020/12/C-Programming-Gate-2020-Set-1-1.jpg"></div>

# 🅲 linguagem C 🔢
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

### Exemplo:

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

### Exemplo:
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

### Tipos de constantes

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

### Exemplo:
```c

"legal"          // constante de string

""               // constante de cadeia nula

"      "             // constante de seis espaços em branco

"A"              // constante de string com caractere único

"Resultado eh\n"     // imprime string com nova linha
```

- Enumerações
A palavra-chave `enum` é usada para definir tipos de enumeração.

### Exemplos:
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
