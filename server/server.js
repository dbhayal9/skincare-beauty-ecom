require('dotenv').config();
const express = require('express');
const cors = require('cors');
const connectDB = require('./config/db');


const app = express();

// Connect Database
connectDB();

// Middleware
app.use(cors());
app.use(express.json());
app.use(express.static('public'));  // Serve static files (images)
// Routes
app.use('/api/products', require('./routes/productRoutes'));

const PORT = process.env.PORT || 5000;
app.listen(PORT, () => console.log(`Server running on port ${PORT}`));
