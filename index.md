# PubMed Article Search Agent with LlamaIndex and OpenAI

## Overview
This project demonstrates an AI-powered agent that efficiently retrieves and summarizes recent PubMed articles based on user queries. Leveraging LlamaIndex, OpenAI, and the Bio.Entrez library, the tool streamlines the process of accessing relevant scientific literature, assisting researchers and professionals in staying up-to-date with cutting-edge findings.

## Features
- **PubMed Integration:** Fetches recent articles related to user queries via the NCBI PubMed database.
- **AI-Powered Summarization:** Uses OpenAI's gpt-3.5-turbo model for intelligent summaries and user interactions.
- **Custom Tools:** Employs LlamaIndex's FunctionTool to create a reusable PubMed querying tool.
- **ReActAgent Framework:** Uses LlamaIndex's ReActAgent for multi-tool querying and contextual understanding.
- **User-Friendly Output:** Displays concise article summaries with titles, authors, publication sources, and dates.

## Architecture
### 1. PubMed Querying
The `query_pubmed` function interacts with the PubMed API to search for articles and retrieve metadata such as:
- Article titles
- Authors
- Publication sources and dates

### 2. Custom Tool Creation
The PubMed querying function is encapsulated into a reusable tool using LlamaIndex's FunctionTool, ensuring seamless integration with the AI agent.

### 3. Agent Initialization
A ReActAgent is instantiated with the PubMed tool and OpenAI's GPT model to handle complex queries and deliver contextual, human-readable responses.

### 4. Interactive Querying
The agent processes natural language queries like:
"Find recent articles on diabetes treatment" and fetches relevant, high-quality articles.

## Technology Stack
- **Programming Language:** Python
- **Libraries:**
  - Bio.Entrez: For PubMed API interactions.
  - LlamaIndex: For creating custom tools and the ReActAgent framework.
  - OpenAI: For GPT-powered natural language interactions.

## How It Works
### 1. Setup
Load the OpenAI API Key and Entrez email securely from environment variables or configuration files.

### 2. Query Execution
The user inputs a natural language query.
The ReActAgent uses the PubMed tool to search for relevant articles and process metadata.

### 3. Response Delivery
The agent summarizes and returns a list of articles, each with:
- Title
- Authors
- Source
- Publication Date

## Example Query and Output
### User Query:
"Find recent articles on diabetes treatment."

### Agent Response:
I found some recent articles related to diabetes treatment:
- **"Comparing Synchronous and Asynchronous Remotely Delivered Lifestyle Interventions: Protocol for a Randomized Noninferiority Trial"**
  - Authors: Pagoto S, Xu R, Bannor R, Idiong C, Goetz J, Fernandes D
  - Source: JMIR Research Protocols
  - Published: December 19, 2024
- **"Multimorbidity in elderly patients with or without T2DM: A real-world cross-sectional analysis based on primary care and hospitalization data"**
  - Authors: Li Y, Geng S, Yuan H, Ge J, Li Q, Chen X, Zhu Y, Liu Y, Guo X, Wang X, Jiang H
  - Source: Journal of Global Health
  - Published: December 20, 2024
- **"Challenges in diagnosis and treatment of KCNJ11-MODY"**
  - Authors: Gonçalves J, Ferreira HU, Ribeiro S, Fernandes da Rocha D, Souto SB, Pedro J, Freitas P, Queirós J
  - Source: Endocrinology, Diabetes & Metabolism Case Reports
  - Published: October 1, 2024

## Installation
### Clone the repository:
```bash
git clone https://github.com/yourusername/PubMed-Search-Agent.git
cd PubMed-Search-Agent
