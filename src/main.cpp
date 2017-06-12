#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}



int main() {
  
  uWS::Hub h;

  // MPC is initialized here
  MPC mpc;

  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"];
	  v *= ( 1609.34 / 3600.0 );
	  
	  
	  double delta = j[1]["steering_angle"];
	  delta *= deg2rad(25);
	  double acceleration = j[1]["throttle"];
	  acceleration *= 1609.34 / 3600.0 / 3600.0;
	  
	  const double Lf = 2.67;

	  // predict state in 100ms plus some extra time for netrwork lag...
	  double latency = 0.11;
	  // so after al this, the car should have a dvanced a little bit (global kinematic model motion of course):
	  double px_lat = px + v * cos(psi)*latency;
	  double py_lat = py + v * sin(psi)*latency;
	  double psi_lat = psi + v*delta/Lf*latency;
	  double v_lat = v + acceleration*latency;
	  
	  // 1. Transforming the control points to the car's local frame (NOTE: WRT original pose reported by the simulator)
	  pair<VectorXd, VectorXd> transformed = TransformPoints2Local(ptsx, ptsy, px, py, psi);
	  VectorXd ptsx_car = transformed.first;
	  VectorXd ptsy_car = transformed.second;
	  
	  // 2. Setup the initial state. NOTE: In local car coordinate frame (original pose reported by the simulator)
	  double x_0 = cos(psi) * ( px_lat - px ) + sin(psi) * (py_lat - py);
	  double y_0 = -sin(psi) * ( px_lat - px ) + cos(psi) * (py_lat - py);
	  double psi_0 = psi_lat - psi;
	  double v_0 = v_lat;
	  // Mow we need cte0 and epsi0. A good chance to fit the polynomial
	  // Ok, now fit a polynomial on the LOCAL points
	  VectorXd coeffs = polyfit(ptsx_car, ptsy_car, 3);
	  
	  // compute the orientation error
	  double cte_0 = polyeval(coeffs, x_0) - y_0;
	  // and the angular error
	  double epsi_0 = psi_0 - atan( polyderiveval(coeffs, x_0) );
	  
	  Eigen::VectorXd state(6);
	  state[0] = x_0; 
	  state[1] = y_0; 
	  state[2] = psi_0; 
	  state[3] = v_0;
	  state[4] = cte_0;
	  state[5] = epsi_0;
	  
	  // solve.......
	  auto vars = mpc.Solve(state, coeffs);
          // store the solver's optimizers
	  double steer_value = vars[vars.size() - 2];
          double throttle_value = vars[vars.size() - 1];
	  
	  // just dump it out for whatever reason...
	  cout <<"New Steering: " << steer_value<<endl;
	  cout <<"New Throttle: " << throttle_value<<endl;
	  
	  
          json msgJson;
          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
          msgJson["steering_angle"] = -steer_value / deg2rad(25);
          msgJson["throttle"] = throttle_value; // less gidder with fixed throttle...

          //Display the MPC predicted trajectory 
          vector<double> mpc_x_vals;
          vector<double> mpc_y_vals;
	  int k = 0;
	  while (k < vars.size() - 2) 
	  {
	    mpc_x_vals.push_back(vars[k]);
	    k++;
	    mpc_y_vals.push_back(vars[k]);
	    k++;
	  }
	  
          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Green line

          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;

          //Display the waypoints/reference line
          vector<double> next_x_vals;
          vector<double> next_y_vals;

	  // filling values 
	  for (int i = 0; i < ptsx_car.size(); i++) 
	  {
	    //next_x_vals.push_back(  cos(psi_0) * ( ptsx_car[i] - x_0 ) +  sin(psi_0) * ( ptsy_car[i] - y_0 ) );
	    //next_y_vals.push_back( -sin(psi_0) * ( ptsx_car[i] - x_0 ) +  cos(psi_0) * ( ptsy_car[i] - y_0 ) );
	      next_x_vals.push_back(ptsx_car[i]);
	      next_y_vals.push_back(ptsy_car[i]);
	      
	    
	  }
          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Yellow line

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;


          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          std::cout << msg << std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          this_thread::sleep_for(chrono::milliseconds(100));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
