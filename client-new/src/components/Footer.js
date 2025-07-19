// import React from 'react';
// import { Container } from 'react-bootstrap';

// const Footer = () => {
//   return (
//     <footer className="py-3 bg-dark text-white">
//       <Container>
//         <div className="text-center">
//           &copy; {new Date().getFullYear()} Skincare & Beauty
//         </div>
//       </Container>
//     </footer>
//   );
// };

// export default Footer;

// import React from 'react';
// import { Container, Row, Col } from 'react-bootstrap';
// import { FaInstagram, FaEnvelope, FaMapMarkerAlt } from 'react-icons/fa';

// const Footer = () => {
//   return (
//     <footer className="bg-dark text-white py-4">
//       <Container>
//         <Row>
//           <Col md={4} className="mb-3">
//             <h5>Glow Beauty</h5>
//             <p className="text-muted">
//               <FaMapMarkerAlt className="me-2" />
//               123 Beauty Street, Cosmopolis
//             </p>
//           </Col>
//           <Col md={4} className="mb-3">
//             <h5>Quick Links</h5>
//             <ul className="list-unstyled">
//               <li><a href="/about" className="text-muted">About Us</a></li>
//               <li><a href="/products" className="text-muted">Our Products</a></li>
//               <li><a href="/contact" className="text-muted">Contact</a></li>
//             </ul>
//           </Col>
//           <Col md={4} className="mb-3">
//             <h5>Connect With Us</h5>
//             <div>
//               <a href="https://instagram.com/glowbeauty" className="text-muted me-3">
//                 <FaInstagram size={24} />
//               </a>
//               <a href="mailto:info@glowbeauty.com" className="text-muted">
//                 <FaEnvelope size={24} />
//               </a>
//             </div>
//             <p className="text-muted mt-2">info@glowbeauty.com</p>
//           </Col>
//         </Row>
//         <Row>
//           <Col className="text-center text-muted">
//             <small>&copy; {new Date().getFullYear()} Glow Beauty. All rights reserved.</small>
//           </Col>
//         </Row>
//       </Container>
//     </footer>
//   );
// };

// export default Footer;

import React from 'react';
import { Container, Row, Col } from 'react-bootstrap';
import { FaInstagram, FaEnvelope } from 'react-icons/fa';

const Footer = () => {
  return (
    <footer className="bg-dark text-white py-4 mt-5">
      <Container>
        <Row>
          <Col md={6}>
            <h5>Glow Beauty</h5>
            <p>123 Beauty Street, Cosmopolis</p>
          </Col>
          <Col md={6} className="text-end">
            <a href="https://instagram.com/glowbeauty" className="text-white me-3">
              <FaInstagram size={24} />
            </a>
            <a href="mailto:info@glowbeauty.com" className="text-white">
              <FaEnvelope size={24} />
            </a>
          </Col>
        </Row>
      </Container>
    </footer>
  );
};

export default Footer;