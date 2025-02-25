#!/usr/bin/env python3

import time
from pymavlink import mavutil

def connect_to_vehicle(connection_string, baud=115200):
    """
    Establish connection to the vehicle
    connection_string: e.g., 'udp:127.0.0.1:14550' or '/dev/ttyUSB0'
    baud: baud rate for serial connections
    """
    print(f"Connecting to vehicle on: {connection_string}")
    vehicle = mavutil.mavlink_connection(connection_string, baud=baud)
    
    # Wait for heartbeat to ensure connection
    vehicle.wait_heartbeat()
    print("Heartbeat received! Connection established.")
    return vehicle

def request_data_stream(vehicle):
    """
    Request attitude data stream from the flight controller
    """
    # Request DATA_STREAM at 10Hz
    vehicle.mav.request_data_stream_send(
        vehicle.target_system,
        vehicle.target_component,
        mavutil.mavlink.MAV_DATA_STREAM_ALL,
        10,  # Rate in Hz
        1    # Start sending
    )
    print("Data stream requested.")

def extract_attitude_data(vehicle):
    """
    Extract and print attitude data from MAVLink messages
    """
    try:
        while True:
            # Receive the next message
            msg = vehicle.recv_match(blocking=True)
            
            if msg is None:
                continue
                
            # Check for ATTITUDE message
            if msg.get_type() == 'ATTITUDE':
                # Convert radians to degrees for readability
                roll = msg.roll * 180.0 / 3.14159
                pitch = msg.pitch * 180.0 / 3.14159
                yaw = msg.yaw * 180.0 / 3.14159
                
                # Print attitude data
                print(f"Roll: {roll:.2f}°, Pitch: {pitch:.2f}°, Yaw: {yaw:.2f}°")
                
            # Small delay to prevent overwhelming the console
            time.sleep(0.1)
            
    except KeyboardInterrupt:
        print("\nStopped by user.")
    except Exception as e:
        print(f"An error occurred: {e}")
    finally:
        vehicle.close()
        print("Connection closed.")

def main():
    # Define your connection string here
    # For serial: '/dev/ttyUSB0' (Linux) or 'COM3' (Windows)
    # For UDP: 'udp:127.0.0.1:14550'
    connection_string = '/dev/ttyUSB0'  # Change this as needed
    baud_rate = 115200  # Common baud rate for iNAV
    
    # Connect to the vehicle
    vehicle = connect_to_vehicle(connection_string, baud_rate)
    
    # Request data stream
    request_data_stream(vehicle)
    
    # Extract and display attitude data
    extract_attitude_data(vehicle)

if __name__ == "__main__":
    main()
