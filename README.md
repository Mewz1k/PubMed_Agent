# **PubMed Article Search Agent with LlamaIndex and OpenAI**

## **Overview**
This project showcases an AI-powered agent designed to retrieve and summarize recent PubMed articles based on user queries. Built using LlamaIndex, OpenAI, and the Bio.Entrez library, this tool demonstrates seamless integration of natural language processing (NLP) with scientific data retrieval to assist researchers and professionals in accessing relevant scientific literature efficiently.

---

## **Features**
- **PubMed Integration**: Fetches recent articles related to user queries using the NCBI PubMed database.
- **AI-Powered Summarization**: Leverages OpenAI's `gpt-3.5-turbo` model for user interactions and streamlined data presentation.
- **Custom Tools**: Utilizes LlamaIndex's `FunctionTool` for creating a PubMed querying tool.
- **ReActAgent Framework**: Employs the ReActAgent framework from LlamaIndex for intelligent multi-tool usage.
- **User-Friendly Output**: Displays articles with details such as title, authors, source, and publication date in a clear and concise format.

---

## **Architecture**
1. **PubMed Querying**: 
   - The `query_pubmed` function interacts with the PubMed API to search for articles and retrieve metadata such as titles, authors, and publication dates.
2. **Custom Tool Creation**: 
   - The PubMed querying function is wrapped into a reusable tool using LlamaIndex's `FunctionTool`.
3. **Agent Initialization**: 
   - A ReActAgent is created using the PubMed tool and OpenAI's language model to handle queries and provide contextual responses.
4. **Interactive Querying**:
   - Users input natural language queries (e.g., "Find recent articles on diabetes treatment") to retrieve relevant articles.

---

## **Technology Stack**
- **Programming Language**: Python
- **Key Libraries**:
  - `Bio.Entrez`: For interacting with PubMed API.
  - `LlamaIndex`: For building the ReActAgent framework and creating tools.
  - `OpenAI`: For GPT-powered natural language interactions.

---

## **How It Works**
1. **Setup**:
   - The API key for OpenAI and email for Entrez are securely loaded from environment variables or files.
2. **Query Execution**:
   - The user enters a query in natural language.
   - The `ReActAgent` uses the PubMed tool to fetch and summarize relevant articles.
3. **Response Delivery**:
   - The agent provides a user-friendly summary of the fetched articles, including key details such as titles, authors, sources, and publication dates.

---

## **Example Output**
**User Query**:  
*"Find recent articles on diabetes treatment."*

**Agent Response**:
I found some recent articles related to diabetes treatment:

"Comparing Synchronous and Asynchronous Remotely Delivered Lifestyle Interventions: Protocol for a Randomized Noninferiority Trial" by Pagoto S, Xu R, Bannor R, Idiong C, Goetz J, Fernandes D, published in JMIR Research Protocols on December 19, 2024.

"Multimorbidity in elderly patients with or without T2DM: A real-world cross-sectional analysis based on primary care and hospitalization data" by Li Y, Geng S, Yuan H, Ge J, Li Q, Chen X, Zhu Y, Liu Y, Guo X, Wang X, Jiang H, published in the Journal of Global Health on December 20, 2024.

"Challenges in diagnosis and treatment of KCNJ11-MODY" by Gonçalves J, Ferreira HU, Ribeiro S, Fernandes da Rocha D, Souto SB, Pedro J, Freitas P, Queirós J, published in Endocrinology, Diabetes & Metabolism Case Reports on October 1, 2024.
