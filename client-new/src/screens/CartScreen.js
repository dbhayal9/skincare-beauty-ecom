import React from 'react';
import { Table, Button, Card, Form } from 'react-bootstrap';
import { Link } from 'react-router-dom';
import { useCart } from '../context/CartContext';

const CartScreen = () => {
  const { cartItems } = useCart();

  return (
    <div className='my-3'>
      <h1>Shopping Cart</h1>
      {cartItems.length === 0 ? (
        <Card>
          <Card.Body>
            Your cart is empty <Link to='/'>Go Back</Link>
          </Card.Body>
        </Card>
      ) : (
        <>
          <Table striped bordered hover responsive className='table-sm'>
            <thead>
              <tr>
                <th>Product</th>
                <th>Price</th>
                <th>Quantity</th>
                <th>Subtotal</th>
              </tr>
            </thead>
            <tbody>
              {cartItems.map((item) => (
                <tr key={item._id}>
                  <td>{item.name}</td>
                  <td>${item.price}</td>
                  <td>
                    <Form.Control
                      as='select'
                      value={item.quantity}
                      onChange={(e) => updateCartHandler(item, e.target.value)}
                    >
                      {[...Array(item.countInStock).keys()].map((x) => (
                        <option key={x + 1} value={x + 1}>
                          {x + 1}
                        </option>
                      ))}
                    </Form.Control>
                  </td>
                  <td>${(item.price * item.quantity).toFixed(2)}</td>
                </tr>
              ))}
            </tbody>
          </Table>
          <div className='d-flex justify-content-end'>
            <h4>Total: ${cartItems.reduce((acc, item) => acc + item.price * item.quantity, 0).toFixed(2)}</h4>
          </div>
          <div className='d-flex justify-content-end mt-3'>
            <Button variant='primary'>Proceed to Checkout</Button>
          </div>
        </>
      )}
    </div>
  );
};

export default CartScreen;