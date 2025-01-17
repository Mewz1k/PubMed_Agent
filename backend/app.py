from flask import Flask, request, jsonify
from flask_sqlalchemy import SQLAlchemy
from Bio import Entrez
import openai
import os

# Initialize Flask app
app = Flask(__name__)

# Database configuration
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///chats.db'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
db = SQLAlchemy(app)

# Database model for storing chats
class Chat(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    project_name = db.Column(db.String(50), nullable=False)
    user_message = db.Column(db.Text, nullable=False)
    agent_response = db.Column(db.Text, nullable=False)

# Load API keys from environment variables
openai.api_key = os.getenv("OPENAI_API_KEY")
Entrez.email = os.getenv("ENTREZ_EMAIL")

# Function to query PubMed
def query_pubmed(query, max_results=5):
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()
        id_list = record.get("IdList", [])

        if not id_list:
            return {"error": "No articles found for the query."}

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
        return {"error": f"Error querying PubMed: {str(e)}"}

# Function to summarize articles using OpenAI
def summarize_articles(articles):
    if not articles or "error" in articles:
        return articles.get("error", "No articles to summarize.")
    
    prompt = "Summarize the following articles:\n"
    for article in articles:
        prompt += f"- {article['Title']} ({article['PubDate']})\n"
    
    try:
        response = openai.Completion.create(
            engine="text-davinci-003",
            prompt=prompt,
            max_tokens=150,
            temperature=0.7,
        )
        return response.choices[0].text.strip()
    except Exception as e:
        return f"Error summarizing articles: {e}"

# API endpoint to handle user queries
@app.route('/query', methods=['POST'])
def query():
    try:
        data = request.json
        query = data.get('query', '').strip()
        project_name = data.get('project_name', 'Default Project').strip()

        if not query:
            return jsonify({"error": "Query cannot be empty."}), 400

        # Query PubMed
        articles = query_pubmed(query)

        # Summarize articles
        response_summary = summarize_articles(articles)

        # Save the chat to the database
        new_chat = Chat(
            project_name=project_name,
            user_message=query,
            agent_response=response_summary
        )
        db.session.add(new_chat)
        db.session.commit()

        return jsonify({"response": response_summary, "articles": articles})

    except Exception as e:
        return jsonify({"error": f"Internal Server Error: {str(e)}"}), 500

# API endpoint to retrieve saved chats
@app.route('/chats', methods=['GET'])
def get_chats():
    try:
        project_name = request.args.get('project', 'Default Project').strip()
        chats = Chat.query.filter_by(project_name=project_name).all()

        if not chats:
            return jsonify({"error": "No chats found for the specified project."}), 404

        return jsonify([
            {"user_message": chat.user_message, "agent_response": chat.agent_response}
            for chat in chats
        ])
    except Exception as e:
        return jsonify({"error": f"Internal Server Error: {str(e)}"}), 500

if __name__ == '__main__':
    # Initialize database and run the app
    with app.app_context():
        db.create_all()  # Create the database if it doesn't exist
    app.run(debug=True)
