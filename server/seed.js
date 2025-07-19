require('dotenv').config();
const mongoose = require('mongoose');
const Product = require('./models/Product');
const connectDB = require('./config/db');

connectDB();

const products = [
  {
    name: "Hydrating Facial Cleanser",
    category: "skincare",
    price: 24.99,
    description: "Gentle cleanser that removes impurities without stripping skin",
    images: ["/images/products/serum1.jpg"], // Make sure this matches your image file
    countInStock: 50,
    featured: true
  },
  {
    name: "Summer Floral Dress",
    category: "dress",
    price: 59.99,
    description: "Lightweight floral dress perfect for summer",
    images: ["/images/products/silk1.jpg"], // Make sure this matches your image file
    countInStock: 20,
    size: "M",
    color: "Multicolor"
  }
];

const seedDB = async () => {
  try {
    await Product.deleteMany();
    await Product.insertMany(products);
    console.log('Database seeded successfully!');
    process.exit();
  } catch (err) {
    console.error('Seeding error:', err);
    process.exit(1);
  }
};

seedDB();
