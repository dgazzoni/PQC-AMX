#!/bin/zsh

CPU=$(sysctl -n machdep.cpu.brand_string | sed -e 's/.*\(M[0-9]\).*/\1/')
COOLDOWN_BEGINNING=60
COOLDOWN_MIDDLE=15

# https://stackoverflow.com/a/226724/523079
cat <<EOF
Requirements and recommendations for running this script:
- Run it from the root of the repository using sudo, which is required to access the CPU cycle counter.
- Run it using the caffeinate command to prevent the system going to sleep during the run: sudo caffeinate ./run_benchmarks.sh
- Try to remove as many sources of variability as possible. We recommend closing down all apps, including menubar ones, and turning off WiFi and Bluetooth. If it is a laptop, we recommend running it connected to AC power, and turning off low power mode.

Would you like to continue?
EOF
select yn in "Yes" "No"; do
    case $yn in
        Yes ) break;;
        No ) exit;;
    esac
done

# https://unix.stackexchange.com/a/389406/493812
if [ "$(id -u)" -ne 0 ]; then
        echo 'This script must be run by root'
        exit 1
fi

if [ -d speed_results_${CPU} ]
then
    # https://stackoverflow.com/a/226724/523079
    echo 'Existing results will be deleted. Are you sure?'
    select yn in "Yes" "No"; do
        case $yn in
            Yes ) rm -rf speed_results_${CPU}; break;;
            No ) exit;;
        esac
    done
fi

mkdir speed_results_${CPU}

rm -rf build
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -G Ninja ..
ninja

echo Waiting for the system to cool down for ${COOLDOWN_BEGINNING} seconds...
sleep ${COOLDOWN_BEGINNING}

SCHEMES=("frodokem" "saber")
PARAMETER_SETS_saber=("lightsaber" "saber" "firesaber")
declare -A PARAMETER_SETS=(
    [frodokem]="640_AES,640_SHAKE,976_AES,976_SHAKE,1344_AES,1344_SHAKE"
    [saber]="lightsaber,saber,firesaber"
)
declare -A IMPLS=(
    [frodokem]="ref"
    [saber]="BHK"
)
declare -A VARIANTS=(
    [frodokem]="ref,opt,neon,opt_amx"
    [saber]="neon,amx_polymul,amx_matmul"
)
declare -A SPEED_SUBROUTINE_PREFIXES=(
    [frodokem]="matmul"
    [saber]="matrixvectormulround"
)

for SCHEME ("$SCHEMES[@]")
do
    for PARAMETER_SET ("${(@s:,:)PARAMETER_SETS[${SCHEME}]}")
    do
        for ALLOC in stack mmap
        do
            for IMPL ("${(@s:,:)IMPLS[${SCHEME}]}")
            do
                for VARIANT ("${(@s:,:)VARIANTS[${SCHEME}]}")
                do
                    for BATCHING in "" 4x
                    do
                        if [ "$SCHEME" = "frodokem" ]
                        then
                            SCHEME_PARAMETER_SET=${SCHEME}_${PARAMETER_SET}
                        else
                            SCHEME_PARAMETER_SET=${PARAMETER_SET}
                        fi

                        if [ "$VARIANT" = "" ]
                        then
                            SPEED_EXEC=speed${BATCHING}_${SCHEME_PARAMETER_SET}_${ALLOC}_${IMPL}
                        else
                            SPEED_EXEC=speed${BATCHING}_${SCHEME_PARAMETER_SET}_${ALLOC}_${IMPL}_${VARIANT}
                        fi

                        if [ -f ${SPEED_EXEC} ]
                        then
                            echo "Benchmarking ${SPEED_EXEC}"
                            
                            # https://stackoverflow.com/questions/77711672/performance-of-cpu-only-code-varies-with-executable-file-name
                            cp ${SPEED_EXEC} speed

                            # Run twice to warm up; e.g. macOS needs to verify the code signature and is slower on the first run
                            ./speed > /dev/null
                            ./speed > /dev/null
                            # Actual run
                            ./speed > "../speed_results_${CPU}/${SCHEME}${BATCHING}:${PARAMETER_SET}:${ALLOC}:${IMPL}:${VARIANT}.txt"

                            echo Waiting for the system to cool down for ${COOLDOWN_MIDDLE} seconds...
                            sleep ${COOLDOWN_MIDDLE}

                            PREFIX=${SPEED_SUBROUTINE_PREFIXES[${SCHEME}]}
                            SPEED_SUBROUTINE_EXEC=speed${BATCHING}_${PREFIX}_${SCHEME_PARAMETER_SET}_${ALLOC}_${IMPL}_${VARIANT}

                            if [ -f ${SPEED_SUBROUTINE_EXEC} ]
                            then
                                echo "Benchmarking ${SPEED_SUBROUTINE_EXEC}"

                                # https://stackoverflow.com/questions/77711672/performance-of-cpu-only-code-varies-with-executable-file-name
                                cp ${SPEED_SUBROUTINE_EXEC} speed

                                # Run twice to warm up; e.g. macOS needs to verify the code signature and is slower on the first run
                                ./speed > /dev/null
                                ./speed > /dev/null
                                # Actual run
                                ./speed > "../speed_results_${CPU}/${PREFIX}${BATCHING}:${PARAMETER_SET}:${ALLOC}:${IMPL}:${VARIANT}.txt"

                                echo Waiting for the system to cool down for ${COOLDOWN_MIDDLE} seconds...
                                sleep ${COOLDOWN_MIDDLE}
                            fi
                        fi
                    done
                done
            done
        done
    done
done

for PARAMETER_SET in 640 976 1344
do
    for IMPL in opt neon opt_amx
    do
        SPEED_EXEC=sample_ct_experiment_frodokem_${PARAMETER_SET}_${IMPL}
        echo "Benchmarking ${SPEED_EXEC}"
        
        # https://stackoverflow.com/questions/77711672/performance-of-cpu-only-code-varies-with-executable-file-name
        cp ${SPEED_EXEC} speed

        # Run twice to warm up; e.g. macOS needs to verify the code signature and is slower on the first run
        ./speed > /dev/null
        ./speed > /dev/null

        # Actual run
        ./speed > "../speed_results_${CPU}/sample_ct:${PARAMETER_SET}:${IMPL}.txt"

        echo Waiting for the system to cool down for ${COOLDOWN_MIDDLE} seconds...
        sleep ${COOLDOWN_MIDDLE}
    done

    SPEED_EXEC=speed_sample_frodokem_${PARAMETER_SET}
    echo "Benchmarking ${SPEED_EXEC}"
    
    # https://stackoverflow.com/questions/77711672/performance-of-cpu-only-code-varies-with-executable-file-name
    cp ${SPEED_EXEC} speed

    # Run twice to warm up; e.g. macOS needs to verify the code signature and is slower on the first run
    ./speed > /dev/null
    ./speed > /dev/null

    # Actual run
    ./speed > "../speed_results_${CPU}/sample:${PARAMETER_SET}.txt"

    echo Waiting for the system to cool down for ${COOLDOWN_MIDDLE} seconds...
    sleep ${COOLDOWN_MIDDLE}
done

cd ..

chown -R $(logname) build/ speed_results_${CPU}/
