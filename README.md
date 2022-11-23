# Accessing Tutorial Notebook
## From Mac
1. Open a terminal window in your local computer and connect to the WSx workstation at Haverford
2. Run “jupyter notebook < password >”
3. Run “ jupyter lab —no-browser —port=8888”
4. Open a new terminal window on your local computer
5. Run “ ssh -N -f -L localhost:8888:localhost:8888 USER@ws9.kinsc-cluster.haverford.edu” from the new window
6. In a browser window go to the address “localhost:8888" and enter the password from step 2

## From Windows
**These steps use the program mobaxterm**
1. Ssh to WS
2. Run “jupyter notebook password”
3. Run “jupyter lab —no-browser —port=8899”
4. Click ’tunneling’ tab one top of Mobaxterm toolbar
5. Hit ’New SSH Tunnel’
6. Set up the display following the example from https://ubccr.freshdesk.com/support/solutions/articles/13000059375-creating-a-ssh-tunnel-using-mobaxterm-on-windows
7. Save the tunnel
8. Open a NEW terminal window on Mobaxterm and enter “ssh -N -f -L localhost:8899:localhost:8899 mstempel@ws9.kinsc-cluster.haverford.edu
9. Open “localhost:8899” in a browser window and enter the password from step 2.
