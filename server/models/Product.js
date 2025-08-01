const mongoose = require('mongoose');

const productSchema = new mongoose.Schema({
  name: { type: String, required: true },
  category: { type: String, required: true, enum: ['skincare', 'dress'] },
  price: { type: Number, required: true },
  description: { type: String, required: true },
  images: { type: [String], required: true },
  countInStock: { type: Number, required: true, default: 0 }
}, { timestamps: true });

module.exports = mongoose.model('Product', productSchema);
