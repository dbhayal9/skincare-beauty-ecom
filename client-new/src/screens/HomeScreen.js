// import React, { useEffect, useState } from 'react';
// import { Row, Col } from 'react-bootstrap';
// import axios from 'axios';
// import Product from '../components/Product';

// const HomeScreen = () => {
//   const [products, setProducts] = useState([]);

//   useEffect(() => {
//     const fetchProducts = async () => {
//       const { data } = await axios.get('http://localhost:5000/api/products');
//       setProducts(data);
//     };
//     fetchProducts();
//   }, []);

//   return (
//     <>
//       <h1>Latest Products</h1>
//       <Row>
//         {products.map((product) => (
//           <Col key={product._id} sm={12} md={6} lg={4} xl={3}>
//             <Product product={product} />
//           </Col>
//         ))}
//       </Row>
//     </>
//   );
// };

// export default HomeScreen;
import React, { useState, useEffect } from 'react';  // Add { useState, useEffect }
import { Row, Col, Carousel } from 'react-bootstrap';
import axios from 'axios';
import Product from '../components/Product';


const HomeScreen = () => {
  const [products, setProducts] = useState([]);

  useEffect(() => {
    const fetchProducts = async () => {
      const { data } = await axios.get('/api/products');
      setProducts(data);
    };
    fetchProducts();
  }, []);

  return (
    <>
      <Carousel className="mb-5">
        <Carousel.Item interval={3000}>
          <img
            className="d-block w-100 banner-img"
            src="/images/banner1.jpg"
            alt="First slide"
          />
          <Carousel.Caption>
            <h3>New Summer Collection</h3>
            <p>Discover our organic skincare line</p>
          </Carousel.Caption>
        </Carousel.Item>
        <Carousel.Item interval={3000}>
          <img
            className="d-block w-100 banner-img"
            src="/images/banner2.jpg"
            alt="Second slide"
          />
          <Carousel.Caption>
            <h3>Special Discount</h3>
            <p>20% off on all beauty dresses</p>
          </Carousel.Caption>
        </Carousel.Item>
      </Carousel>

      <h2 className="mb-4">Featured Products</h2>
      <Row>
        {products.map((product) => (
          <Col key={product._id} sm={12} md={6} lg={4} xl={3}>
            <Product product={product} />
          </Col>
        ))}
      </Row>
    </>
  );
};

export default HomeScreen;