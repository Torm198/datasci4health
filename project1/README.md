*2026.1 Ciência e Visualização de Dados em Saúde*

# Projeto: Atlas da Resistência: Análise Comparativa de Redes de Coexpressão Transcriptômica em Patógenos ESKAPE sob Estresse por Carbapenêmicos

# Project: The Resistance Atlas: Comparative Analysis of Transcriptomic Co-expression Networks in ESKAPE Pathogens under Carbapenem Stress

# Descrição Resumida do Projeto

O projeto Atlas da Resistência investiga a resposta adaptativa de patógenos do grupo ESKAPE, especificamente Klebsiella pneumoniae, Acinetobacter baumannii e Pseudomonas aeruginosa, quando submetidos ao estresse pelo antibiótico Meropenem. A motivação central reside no fato de que a resistência antimicrobiana não é um evento isolado de um único gene, mas uma propriedade emergente de sistemas biológicos complexos que se organizam para garantir a sobrevivência bacteriana. No contexto clínico atual, essas três bactérias representam as ameaças mais críticas em ambientes hospitalares devido à sua capacidade de "escapar" da ação de carbapenêmicos, que são frequentemente a última linha de defesa terapêutica. O problema abordado pelo projeto é a falta de uma visão sistêmica e comparativa que identifique se diferentes espécies utilizam uma arquitetura de rede comum para resistir ao mesmo fármaco. Utilizando dados de transcriptoma (RNA-Seq) obtidos de bases públicas, o trabalho emprega a Ciência de Redes para transformar níveis de expressão gênica em grafos de coexpressão, onde os genes atuam como nós e suas correlações funcionais como arestas. O objetivo final é realizar uma análise visual e topológica para identificar genes-hub e módulos de resistência conservados, permitindo determinar se existe um "core" transcriptômico universal que possa ser explorado como alvo para novas estratégias de tratamento que ignorem as fronteiras entre espécies.

# Slides

> Coloque aqui o link para o PDF da apresentação da parte 1.

# Fundamentação Teórica

A base biológica deste projeto reside no estresse da parede celular induzido pelos carbapenêmicos, fármacos que inativam as proteínas ligadoras de penicilina (PBPs), impedindo a síntese de peptidoglicano e levando à lise bacteriana (ou seja, o antibiótico impede a bactéria de manter sua “parede protetora”, fazendo com que ela se rompa e morra).
Em resposta, patógenos Gram-negativos ativam uma cascata transcricional complexa (uma série de reações coordenadas no funcionamento interno da célula) que envolve a repressão de porinas (redução de “portas” por onde o antibiótico entra) para limitar a entrada da droga, a superexpressão de bombas de efluxo (ativação de mecanismos que “expulsam” o antibiótico para fora) para sua remoção e a reorganização do metabolismo energético (ajustes no uso de energia da célula) para sustentar mecanismos de reparo (processos que tentam consertar os danos causados).
A análise de redes é fundamental neste contexto, pois a resistência não depende apenas da presença de genes como blaKPC ou blaNDM (genes conhecidos por conferir resistência), mas da interação dinâmica entre reguladores globais e vias metabólicas secundárias (diferentes partes do sistema celular que se comunicam e se influenciam) que, em conjunto, conferem o fenótipo de multirresistência (ou seja, a capacidade da bactéria de resistir a vários antibióticos ao mesmo tempo).

# Perguntas de Pesquisa

1. Quais genes aumentam sua expressão simultaneamente nas três bactérias quando elas tentam sobreviver ao Meropenem?
2. Se transformarmos esses dados em um grafo, quais são os 5 genes mais conectados (hubs) que aparecem como "líderes" da resistência em cada espécie?
3. isualmente, as redes de resistência dessas três bactérias se parecem ou cada uma tem um "estilo" de defesa completamente diferente?


# Bases de Dados

| Base de Dados | Endereço na Web | Resumo descritivo |
|-|-|-|
| NCBI Gene Expression Omnibus (GEO)                        | https://www.ncbi.nlm.nih.gov/geo/ | O maior repositório público de dados de genômica funcional do mundo. Permite buscar especificamente por experimentos de RNA-Seq de patógenos sob estresse de antibióticos.           |
| BV-BRC (Bacterial & Viral Bioinformatics Resource Center) | https://www.bv-brc.org/           | Antigo PATRIC, é uma base especializada em patógenos bacterianos. Oferece ferramentas integradas para comparar transcriptomas de diferentes cepas de Klebsiella e Pseudomonas.       |
| European Nucleotide Archive (ENA)                         | https://www.ebi.ac.uk/ena/        | Repositório europeu que armazena sequências brutas e processadas. É uma alternativa essencial caso algum dataset relevante de parceiros internacionais não esteja espelhado no NCBI. |
| PubMed Central (PMC)                                      | https://www.ncbi.nlm.nih.gov/pmc/ | Base de artigos científicos de acesso aberto. Essencial para baixar as "tabelas suplementares" de artigos recentes que contêm os valores de Fold-Change já calculados pelos autores. |
| The Comprehensive Antibiotic Resistance Database          | https://card.mcmaster.ca/         | Base de referência para identificar genes de resistência. Será usada para anotar os "nós" do grafo e confirmar se os hubs encontrados possuem função conhecida de resistência.       |

# Modelo Lógico

O modelo lógico descreve como os dados do projeto são estruturados e relacionados para permitir a construção e análise de redes de coexpressão gênica.

## 1. Entidades

### Expressão Gênica
Representa cada gene analisado no estudo.

**Atributos:**
- `gene_id` (PK): Identificador único do gene
- `nome_gene`: Nome do gene
- `organismo`: Bactéria de origem
- `anotacao_resistencia`: Indica se o gene está associado à resistência (ex: presente no CARD)
- `funcao_biologica`: Função biológica do gene
- `amostra_id`: Identificador da amostra associada
- `nivel_expresao`: Nível de expressão do gene

---

### Condição Experimental
Define o contexto em que as amostras foram coletados.

**Atributos:**
- `amostra_id`: Identificador da amostra
- `antibiotico`: Condição experimental (ex: meropenem ou normal)

---

### Coexpressão
Define a relação associada à expressão gênica em diferentes condições.

**Atributos:**
- `gene_id`: Identificador do gene
- `nivel_expresao_meropenem`: Nível de expressão sob presença de antibiótico
- `nivel_expresao_normal`: Nível de expressão em condição normal
- `diferenca_expressao`: Diferença entre os níveis de expressão

---

### Interações entre Genes (STRING)
Representa as conexões entre pares de genes/proteínas com base na base STRING, incluindo interações conhecidas e preditas.
Atributos:

**Atributos:**
- `gene_id_1`: Identificador do primeiro gene
- `gene_id_2`: Identificador do segundo gene
- `score_interacao`: Score combinado de confiança da interação (0 a 1)
- `tipo_interacao`: Tipo da interação (ex: experimental, coexpressão, banco de dados, text mining)
- `evidencia_experimental`: Score baseado em experimentos laboratoriais
- `evidencia_coexpressao`: Score baseado em padrões de expressão similares
- `evidencia_textmining`: Score baseado em coocorrência em artigos científicos
- `evidencia_database`: Score baseado em bases curadas

---

## 2. Aplicação às Perguntas do Projeto

- **Genes coexpressos entre espécies:** análise de **Expressão Gênica** + **Coexpressão**
- **Identificação de hubs:** uso de **Interações entre Genes** + **Coexpressão**

---

# Metodologia
A execução do projeto será dividida em quatro etapas principais, integrando a análise biológica de expressão gênica com a modelagem matemática de redes:

## 1. Coleta e Pré-processamento de Dados

Utilizaremos dados brutos de contagem (*read counts*) de sequenciamento de RNA (RNA-Seq) provenientes do repositório público NCBI Gene Expression Omnibus (GEO). Serão selecionados três datasets independentes para:

- *Klebsiella pneumoniae*  
- *Acinetobacter baumannii*  
- *Pseudomonas aeruginosa*  

Todos contendo:
- grupos controle (sem tratamento)  
- grupos submetidos ao estresse por Meropenem  

---

## 2. Análise de Expressão Diferencial (Fold-Change)

Para identificar os genes que respondem ao antibiótico, utilizaremos pacotes estatísticos como:

- DESeq2  
- EdgeR  

O critério de seleção será baseado no **log₂ Fold-Change (LFC)**:

- Filtrar genes com:  
  - |LFC| ≥ 2.0  
  - P-valor < 0.05  

Esse filtro garante que a rede represente uma resposta biológica real ao estresse, reduzindo ruídos metabólicos basais.

---

## 3. Construção da Rede de Coexpressão

A partir dos genes filtrados, serão construídas matrizes de adjacência baseadas em:

- Correlação de Pearson  
- Correlação de Spearman  

**Definições da rede:**

- **Nós:** genes diferencialmente expressos  
- **Arestas:** conexões entre genes com correlação acima de um limiar crítico (*threshold*), indicando co-regulação sob estresse  

---

## 4. Análise de Redes e Visualização

Serão aplicadas técnicas de análise de redes, incluindo:

- **Centralidade de Grau:** identificação de genes-hub (nós mais conectados)  
- **Detecção de Comunidades (Louvain):** agrupamento de genes em módulos funcionais  

A visualização será realizada em ferramentas como:

- Cytoscape  
- Gephi  

**Configuração visual:**
- Cor dos nós:  
  - Vermelho → aumento de expressão  
  - Azul → diminuição de expressão  
- Tamanho dos nós: proporcional à centralidade (importância na rede)  

---

# Ferramentas

> Ferramentas a serem utilizadas (com base na visão atual do grupo sobre o projeto).

# Referências Bibliográficas

## Referências Bibliográficas

CASTANHEIRA, Mariana et al. The plethora of resistance mechanisms in *Pseudomonas aeruginosa*: transcriptome analysis reveals a potential role of lipopolysaccharide pathway proteins to novel β-lactam/β-lactamase inhibitor combinations.

ABBAS, Fatima Moeen; AJEEL, Murtadha Hamza. Carbapenemases, types and epidemiology: a review. 2023.

GEISINGER, Edward; ISBERG, Ralph R. Interplay between antibiotic resistance and virulence during disease promoted by multidrug-resistant bacteria.

MARTÍNEZ-MARTÍNEZ, Luis. Carbapenemases: the never-ending story.

QIN, Hao et al. Comparative transcriptomics of multidrug-resistant *Acinetobacter baumannii* in response to antibiotic treatments.