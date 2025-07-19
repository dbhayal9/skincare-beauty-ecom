// import React from 'react';
// import { Container } from 'react-bootstrap';
// import { BrowserRouter as Router, Route, Routes } from 'react-router-dom';
// import Header from './components/Header';
// import Footer from './components/Footer';
// import HomeScreen from './screens/HomeScreen';
// import ProductScreen from './screens/ProductScreen';

// // Add these routes




// function App() {
//   return (
//     <Router>
//       <Header />
//       <main className="py-3">
//         <Container>
//           <Routes>
//             <Route path="/" element={<HomeScreen />} exact />
//             <Route path="/product/:id" element={<ProductScreen />} />
//             <Route path="/about" element={<AboutScreen />} />
//             <Route path="/contact" element={<ContactScreen />} />
//             <Route path="/products" element={<ProductListScreen />} />
//           </Routes>
//         </Container>
//       </main>
//       <Footer />
//     </Router>
//   );
// }

// // export default App;
// import React from 'react';
// import { Container } from 'react-bootstrap';
// import { BrowserRouter as Router, Routes, Route } from 'react-router-dom';
// import Header from './components/Header';
// import Footer from './components/Footer';
// import HomeScreen from './screens/HomeScreen';
// import ProductScreen from './screens/ProductScreen';
// import AboutScreen from './screens/AboutScreen';
// import ContactScreen from './screens/ContactScreen';


// function App() {
//   return (
//     <Router>
//       <Header />
//       <main className="py-3">
//         <Container>
//           <Routes>
//             <Route path="/" element={<HomeScreen />} exact />
//             <Route path="/products" element={<HomeScreen />} />
//             <Route path="/about" element={<AboutScreen />} />
//             <Route path="/contact" element={<ContactScreen />} />
//             <Route path="/product/:id" element={<ProductScreen />} />
//           </Routes>
//         </Container>
//       </main>
//       <Footer />
//     </Router>
//   );
// }

// export default App;

// export default App;
import React from 'react';
import { Container } from 'react-bootstrap';
import { BrowserRouter as Router, Routes, Route } from 'react-router-dom';
import Header from './components/Header';
import Footer from './components/Footer';
import HomeScreen from './screens/HomeScreen';
import ProductScreen from './screens/ProductScreen';
import AboutScreen from './screens/AboutScreen';
import ContactScreen from './screens/ContactScreen';
import { CartProvider } from './context/CartContext';
import CartScreen from './screens/CartScreen';

function App() {
  return (
    <CartProvider>
      <Router>
        <Header />
        <main className="py-3">
          <Container>
            <Routes>
              <Route path="/" element={<HomeScreen />} exact />
              <Route path="/products" element={<HomeScreen />} />
              <Route path="/about" element={<AboutScreen />} />
              <Route path="/contact" element={<ContactScreen />} />
              <Route path="/product/:id" element={<ProductScreen />} />
              <Route path='/cart' element={<CartScreen />} />
            </Routes>
          </Container>
        </main>
        <Footer />
      </Router>
    </CartProvider>
  );
}

export default App;