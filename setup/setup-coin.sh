#!/bin/bash


# Function to display a message
log_message() {
    echo "===================================================================="
    echo "$1"
    echo "===================================================================="
}

# Update system packages
#log_message "Updating system packages..."
#sudo apt-get update -y
#sudo apt-get install -y build-essential clang git wget

# Installation directory
DEFAULT_INSTALL_DIR="/opt"
DIR_NAME="/coin-or"
echo "The default installation directory for the needed COIN-OR Tools is: $DEFAULT_INSTALL_DIR$DIR_NAME"
echo "Press Enter to accept this default location, or specify a different directory:"
read -r -e -p ":" USER_INPUT
if [ -n "$USER_INPUT" ]; then
    INSTALL_DIR="$USER_INPUT"
else
    INSTALL_DIR=${DEFAULT_INSTALL_DIR}
fi

echo "The required COIN-OR tools will be installed in: $INSTALL_DIR"
export EXTERN_SOLVER_DIR=${INSTALL_DIR}$DIR_NAME

if [ -d "${INSTALL_DIR}${DIR_NAME}" ]; then
    echo "There already exists a folder called coin-or in ${INSTALL_DIR}. Skip installation."
    return
fi

run_cmd() {
    if [ -w "${INSTALL_DIR}${DIR_NAME}" ]; then
        "$@"
    else
        sudo "$@"
    fi
}

if [ ! -d "${INSTALL_DIR}${DIR_NAME}" ]; then
  echo "The directory ${INSTALL_DIR}${DIR_NAME} does not exist. Creating it now..."
  # Try creating the directory
  mkdir -vp "${INSTALL_DIR}${DIR_NAME}"  # Without sudo first
  # Check if the directory creation failed due to permission issues
  if [ $? -ne 0 ]; then
    # If permission error, then use sudo
    echo "Error: Failed to create directory ${INSTALL_DIR}${DIR_NAME}. Trying with sudo..."
    run_cmd mkdir -vp "${INSTALL_DIR}${DIR_NAME}" || {
      echo "Error: Failed to create directory ${INSTALL_DIR}${DIR_NAME} even with sudo."
      exit 1
    }
  fi
fi

if [ ! -w "${INSTALL_DIR}${DIR_NAME}" ]; then
    echo "You need sudo to write to ${INSTALL_DIR}${DIR_NAME}."
fi
cd ${INSTALL_DIR}${DIR_NAME}
run_cmd wget https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew
run_cmd chmod u+x coinbrew
run_cmd ./coinbrew fetch Cbc
run_cmd ./coinbrew build Cbc --tests none

log_message "All COIN-OR libraries installed successfully in ${INSTALL_DIR}${DIR_NAME}"
