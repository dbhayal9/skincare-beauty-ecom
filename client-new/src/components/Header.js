// import React from 'react';
// import { Navbar, Container } from 'react-bootstrap';

// const Header = () => {
//   return (
//     <Navbar bg="dark" variant="dark">
//       <Container>
//         <Navbar.Brand href="/">Skincare & Beauty</Navbar.Brand>
//       </Container>
//     </Navbar>
//   );
// };

// export default Header;

// import React from 'react';
// import { Navbar, Nav, Container } from 'react-bootstrap';
// import { Link } from 'react-router-dom';

// const Header = () => {
//   return (
//     <Navbar bg="dark" variant="dark" expand="lg" collapseOnSelect>
//       <Container>
//         <Link to="/" className="navbar-brand">
//           Glow Beauty
//         </Link>
//         <Navbar.Toggle aria-controls="basic-navbar-nav" />
//         <Navbar.Collapse id="basic-navbar-nav">
//           <Nav className="ms-auto">
//             <Link to="/" className="nav-link">
//               Home
//             </Link>
//             <Link to="/products" className="nav-link">
//               Products
//             </Link>
//             <Link to="/about" className="nav-link">
//               About
//             </Link>
//             <Link to="/contact" className="nav-link">
//               Contact
//             </Link>
//           </Nav>
//         </Navbar.Collapse>
//       </Container>
//     </Navbar>
//   );
// };

// export default Header;
import React from 'react';
import { Navbar, Nav, Container } from 'react-bootstrap';
import { Link } from 'react-router-dom';
import { FaShoppingCart } from 'react-icons/fa';

const Header = () => {
  return (
    <Navbar bg="dark" variant="dark" expand="lg">
      <Container>
        <Link to="/" className="navbar-brand">Sheenel</Link>
        <Navbar.Toggle aria-controls="basic-navbar-nav" />
        <Navbar.Collapse id="basic-navbar-nav">
          <Nav className="ms-auto">
            <Link to="/" className="nav-link">Home</Link>
            <Link to="/products" className="nav-link">Products</Link>
            <Link to="/about" className="nav-link">About</Link>
            <Link to="/contact" className="nav-link">Contact</Link>
          </Nav>
        </Navbar.Collapse>
      </Container>
    </Navbar>
  );
};

export default Header;