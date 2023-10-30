//
// SPDX-License-Identifier: LGPL-3.0-or-later
//
import inet.networklayer.configurator.ipv4.Ipv4NetworkConfigurator;
import inet.node.contract.INetworkNode;
import inet.physicallayer.wireless.common.contract.packetlevel.IRadioMedium;
import inet.visualizer.contract.IIntegratedVisualizer;
import inet.networklayer.configurator.ipv4.Ipv4NetworkConfigurator;
import inet.node.ethernet.Eth100M;
import inet.node.inet.SensorNode;
import inet.node.inet.StandardHost;
import inet.physicallayer.wireless.apsk.packetlevel.ApskScalarRadioMedium;
import inet.visualizer.common.IntegratedVisualizer;

network topology
{
    @display("bgb=100,100");
    submodules:
        configurator: Ipv4NetworkConfigurator {
            @display("p=100,100");
        }
        visualizer: IntegratedVisualizer {
            @display("p=4.0722203,9.346079");
        }
        radioMedium: ApskScalarRadioMedium {
            @display("p=4.0054626,14.219392");
        }
        
        sensor0:SensorNode {
           @display("p=74.2,79.9;i=misc/sensor2");
        }
        sensor1:SensorNode {
           @display("p=89.2,53.5;i=misc/sensor2");
        }
        sensor2:SensorNode {
           @display("p=11.4,89.5;i=misc/sensor2");
        }
        sensor3:SensorNode {
           @display("p=73.0,41.6;i=misc/sensor2");
        }
        sensor4:SensorNode {
           @display("p=71.6,42.1;i=misc/sensor2");
        }
        sensor5:SensorNode {
           @display("p=59.2,18.8;i=misc/sensor2");
        }
        sensor6:SensorNode {
           @display("p=55.6,47.4;i=misc/sensor2");
        }
        sensor7:SensorNode {
           @display("p=21.3,82.2;i=misc/sensor2");
        }
        sensor8:SensorNode {
           @display("p=71.1,87.6;i=misc/sensor2");
        }
        sensor9:SensorNode {
           @display("p=19.2,98.4;i=misc/sensor2");
        }
        sensor10:SensorNode {
           @display("p=5.8,23.3;i=misc/sensor2");
        }
        sensor11:SensorNode {
           @display("p=26.8,46.1;i=misc/sensor2");
        }
        sensor12:SensorNode {
           @display("p=52.0,91.8;i=misc/sensor2");
        }
        sensor13:SensorNode {
           @display("p=44.4,70.4;i=misc/sensor2");
        }
        sensor14:SensorNode {
           @display("p=20.3,8.6;i=misc/sensor2");
        }
        sensor15:SensorNode {
           @display("p=39.9,90.9;i=misc/sensor2");
        }
        sensor16:SensorNode {
           @display("p=96.7,74.9;i=misc/sensor2");
        }
        sensor17:SensorNode {
           @display("p=75.2,70.2;i=misc/sensor2");
        }
        sensor18:SensorNode {
           @display("p=43.1,67.4;i=misc/sensor2");
        }
        sensor19:SensorNode {
           @display("p=5.2,14.2;i=misc/sensor2");
        }
        sensor20:SensorNode {
           @display("p=21.0,79.3;i=misc/sensor2");
        }
        sensor21:SensorNode {
           @display("p=89.1,78.3;i=misc/sensor2");
        }
        sensor22:SensorNode {
           @display("p=17.3,78.0;i=misc/sensor2");
        }
        sensor23:SensorNode {
           @display("p=30.4,41.9;i=misc/sensor2");
        }
        sensor24:SensorNode {
           @display("p=60.6,42.3;i=misc/sensor2");
        }
        sensor25:SensorNode {
           @display("p=66.5,35.9;i=misc/sensor2");
        }
        sensor26:SensorNode {
           @display("p=94.0,76.0;i=misc/sensor2");
        }
        sensor27:SensorNode {
           @display("p=22.8,43.2;i=misc/sensor2");
        }
        sensor28:SensorNode {
           @display("p=64.4,38.2;i=misc/sensor2");
        }
        sensor29:SensorNode {
           @display("p=26.7,16.7;i=misc/sensor2");
        }
        
        sensor30:SensorNode {
           @display("p=22.8,18.6;i=misc/sensor2");
        }
        sensor31:SensorNode {
           @display("p=99.9,19.6;i=misc/sensor2");
        }
        
        
 
        gateway: SensorNode {
            @display("p=23.89926,17.557278;i=misc/sensorgateway");
        }
        server: StandardHost {
            @display("p=35.381584,17.557278");
        }
    connections:
        gateway.ethg++ <--> Eth100M <--> server.ethg++;
}

