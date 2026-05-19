#

## Commands

```bash
# Initiates project.
cmake -B build -S .

# Removes current license path.
unset GRB_LICENSE_FILE
# Creates system path to the license.
export GRB_LICENSE_FILE="example_path"

# Edit file for chaning license path.
nano ~/.zshrc

# Builds app.
cmake --build build

# Runs app.
./build/EVCSP_Application
```