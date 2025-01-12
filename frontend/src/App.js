import React, { useState } from 'react';
import './App.css';

function App() {
  const [query, setQuery] = useState('');
  const [chats, setChats] = useState([]);

  const handleSubmit = async (e) => {
    e.preventDefault();
    const userMessage = { role: 'user', content: query };
    setChats([...chats, userMessage]);

    const response = await fetch('http://localhost:5000/query', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ query }),
    });

    const data = await response.json();
    const agentMessage = { role: 'agent', content: data.response };
    setChats([...chats, userMessage, agentMessage]);
    setQuery('');
  };

  return (
    <div className="chat-container">
      <div className="chat-box">
        {chats.map((chat, index) => (
          <div key={index} className={chat.role}>
            {chat.content}
          </div>
        ))}
      </div>
      <form onSubmit={handleSubmit} className="input-container">
        <input
          type="text"
          value={query}
          onChange={(e) => setQuery(e.target.value)}
          placeholder="Ask me anything..."
        />
        <button type="submit">Send</button>
      </form>
    </div>
  );
}

export default App;
