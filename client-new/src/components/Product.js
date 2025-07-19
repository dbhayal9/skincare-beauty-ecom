// import React from 'react';
// import { Card } from 'react-bootstrap';
// import { Link } from 'react-router-dom';

// const Product = ({ product }) => {
//   return (
//     <Card className='my-3 p-3 rounded'>
//       <Link to={`/product/${product._id}`}>
//         <Card.Img
//           src={`http://localhost:5000${product.images[0]}`}
//           variant="top"
//           onError={(e) => {
//             e.target.src = 'https://via.placeholder.com/300';
//             console.error('Image load failed:', product.images[0]);
//           }}
//         />
//       </Link>

//       <Card.Body>
//         <Link to={`/product/${product._id}`}>
//           <Card.Title as='div'>
//             <strong>{product.name}</strong>
//           </Card.Title>
//         </Link>
//         <Card.Text as='h3'>${product.price}</Card.Text>
//       </Card.Body>
//     </Card>
//   );
// };

// export default Product;

// worked
// import React, { useState } from 'react';
// import { Card, Carousel } from 'react-bootstrap';
// import { Link } from 'react-router-dom';

// const Product = ({ product }) => {
//   const [index, setIndex] = useState(0);

//   return (
//     <Card className='my-3 p-3 rounded product-card'>
//       <Carousel 
//         activeIndex={index} 
//         onSelect={setIndex}
//         interval={null}
//         indicators={product.images.length > 1}
//       >
//         {product.images.map((img, i) => (
//           <Carousel.Item key={i}>
//             <Link to={`/product/${product._id}`}>
//               <img
//                 src={`http://localhost:5000${img}`}
//                 className="d-block w-100 product-image"
//                 alt={`${product.name} angle ${i+1}`}
//               />
//             </Link>
//           </Carousel.Item>
//         ))}
//       </Carousel>
      
//       <Card.Body>
//         <Link to={`/product/${product._id}`}>
//           <Card.Title as='div' className='product-title'>
//             <strong>{product.name}</strong>
//           </Card.Title>
//         </Link>
//         <Card.Text as='h3' className='product-price'>${product.price}</Card.Text>
//       </Card.Body>
//     </Card>
//   );
// };

// export default Product;

import React, { useState } from 'react';
import { Card, Form, Button } from 'react-bootstrap';
import { Link } from 'react-router-dom';
import { useCart } from '../context/CartContext';
const Product = ({ product }) => {
  const { addToCart } = useCart();
  const [quantity, setQuantity] = useState(1);

  const addToCartHandler = () => {
    addToCart(product, quantity);
  };

  return (
    <Card className='my-3 p-3 rounded'>
      <Link to={`/product/${product._id}`}>
        <Card.Img src={`http://localhost:5000${product.images[0]}`} variant='top' />
      </Link>

      <Card.Body>
        <Link to={`/product/${product._id}`}>
          <Card.Title as='div'>
            <strong>{product.name}</strong>
          </Card.Title>
        </Link>

        <Card.Text as='h3'>${product.price}</Card.Text>

        {product.countInStock > 0 && (
          <>
            <Form.Group controlId='quantity' className='my-2'>
              <Form.Label>Quantity</Form.Label>
              <Form.Control
                as='select'
                value={quantity}
                onChange={(e) => setQuantity(Number(e.target.value))}
              >
                {[...Array(product.countInStock).keys()].map((x) => (
                  <option key={x + 1} value={x + 1}>
                    {x + 1}
                  </option>
                ))}
              </Form.Control>
            </Form.Group>

            <Button
              onClick={addToCartHandler}
              className='btn-block w-100'
              disabled={product.countInStock === 0}
            >
              {product.countInStock > 0 ? 'Add to Cart' : 'Out of Stock'}
            </Button>
          </>
        )}
      </Card.Body>
    </Card>
  );
};

export default Product;