from flask import Flask, request, jsonify
from flask_sqlalchemy import SQLAlchemy
from Bio import Entrez
import openai

app = Flask(Research Agent)
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///chats.db'
db = SQLAlchemy(app)

# Database model for storing chats
class Chat(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    project_name = db.Column(db.String(50), nullable=False)
    user_message = db.Column(db.Text, nullable=False)
    agent_response = db.Column(db.Text, nullable=False)

# Load API keys (upload later)
openai.api_key = "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"  # Replace with your API key
Entrez.email = "chayyalavarthi@gmail.com"  # Replace with your email

# Function to query PubMed
def query_pubmed(query, max_results=5):
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()
        id_list = record["IdList"]

        articles = []
        for pubmed_id in id_list:
            summary_handle = Entrez.esummary(db="pubmed", id=pubmed_id, retmode="xml")
            summary_record = Entrez.read(summary_handle)
            summary_handle.close()
            article = summary_record[0]
            articles.append({
                "Title": article.get('Title', 'N/A'),
                "Authors": ', '.join(article.get('AuthorList', [])),
                "Source": article.get('Source', 'N/A'),
                "PubDate": article.get('PubDate', 'N/A')
            })
        return articles
    except Exception as e:
        return {"error": str(e)}

# API endpoint to handle user queries
@app.route('/query', methods=['POST'])
def query():
    data = request.json
    query = data.get('query', '')
    project_name = data.get('project_name', 'Default Project')
    
    # Query PubMed
    articles = query_pubmed(query)

    # Use OpenAI GPT to summarize articles (mock example here)
    response_summary = f"Summarized response for query: {query}"
    
    # Save the chat to the database
    new_chat = Chat(project_name=project_name, user_message=query, agent_response=response_summary)
    db.session.add(new_chat)
    db.session.commit()

    return jsonify({"response": response_summary, "articles": articles})

# API endpoint to retrieve saved chats
@app.route('/chats', methods=['GET'])
def get_chats():
    project_name = request.args.get('project', 'Default Project')
    chats = Chat.query.filter_by(project_name=project_name).all()
    return jsonify([{"user_message": chat.user_message, "agent_response": chat.agent_response} for chat in chats])

if __name__ == '__main__':
    db.create_all()  # Create the database if it doesn't exist
    app.run(debug=True)
