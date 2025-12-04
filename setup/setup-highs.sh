DEFAULT_INSTALL_DIR="/opt"
DIR_NAME="/HiGHS"

echo "The default installation directory for HiGHS is: $DEFAULT_INSTALL_DIR/HiGHS"
echo "Press Enter to accept this default location, or specify a different directory:"
read -r -e -p ":" USER_INPUT
if [ -n "$USER_INPUT" ]; then
    INSTALL_DIR="$USER_INPUT"
else
    INSTALL_DIR=${DEFAULT_INSTALL_DIR}
fi

echo "HiGHS will be installed in: $INSTALL_DIR"
export EXTERN_SOLVER_DIR=${INSTALL_DIR}$DIR_NAME

if [ -d "${INSTALL_DIR}${DIR_NAME}" ]; then
    echo "There already exists a HiGHS folder in ${INSTALL_DIR}. Skip installation."
    return
fi

run_cmd() {
    if [ -w "$INSTALL_DIR" ]; then
        "$@"
    else
        sudo "$@"
    fi
}

if [ ! -d "$INSTALL_DIR" ]; then
  echo "The directory $INSTALL_DIR does not exist. Creating it now..."
  # Try creating the directory
  mkdir -vp "$INSTALL_DIR"  # Without sudo first
  # Check if the directory creation failed due to permission issues
  if [ $? -ne 0 ]; then
    # If permission error, then use sudo
    echo "Error: Failed to create directory ${INSTALL_DIR}. Trying with sudo..."
    run_cmd mkdir -vp "$INSTALL_DIR" || {
      echo "Error: Failed to create directory ${INSTALL_DIR} even with sudo."
      exit 1
    }
  fi
fi
if [ ! -w "$INSTALL_DIR" ]; then
    echo "You need sudo to write to $INSTALL_DIR."
fi
cd $INSTALL_DIR
run_cmd git clone https://github.com/ERGO-Code/HiGHS.git
cd HiGHS
run_cmd cmake -B build
run_cmd cmake --build build

