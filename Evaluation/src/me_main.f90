! ****************************************
! Main program for emulator training and testing
! Modified 11.2.2025 by Tomi Raatikainen (FMI)
! ****************************************

PROGRAM me_main
    ! Emulator training and testing modules
    USE me_training, ONLY : leave_one_out, leave_some_out, test_emulator
    ! Tests with IFS data
    USE me_ifs_sample, ONLY : IFS_sample_full
    USE me_ifs_day, ONLY : IFS_test_day
    !
    ! Emulator development and tests
    ! ======================
    ! Train emulator first (and save emulator) and then test by predicting updraft velocity for each input
    ! Fixed random number generator seed allows the production of the same emulator 
    !
    ! Tests with IFS data
    !-----------------------
    ! Determine emulator inputs and predict updafts for the sample dataset
    !CALL IFS_sample_full(output='test_emu_sample_all_N24.dat')
    !CALL IFS_sample_full(output='test_emu_sample_all_J25.dat')
    CALL IFS_sample_full(output='test_emu_sample_all_Sc_J26.dat')
    !CALL IFS_sample_full(output='test_emu_sample_all_H26.dat')
    !
    ! Complete tests
    ! -----------------
    IF (.FALSE.) THEN
        ! Data from January 2025 - compile with "j25" flag
        ! 1) Training
        CALL test_emulator('DATA_wp.dat',' ','test_wp_J25.dat','gp_wp_J25.emu',iseed=12)
        ! 2)  Leave-n-out, where n=0.1*ntot
        CALL leave_some_out(0.1,'DATA_wp.dat','gp_wp_J25.emu','lnu_wp_J25.dat',iseed=12)
        ! 3)  IFS predictions
        CALL IFS_sample_full(output='test_emu_sample_all_J25.dat')
    ELSEIF (.FALSE.) THEN
        ! Data from January 2026 (Stratocumulus only) - compile with "j26" flag
        ! 1) Training
        CALL test_emulator('DATA_ws_Sc.dat',' ','test_ws_Sc_J26.dat','gp_ws_Sc_J26.emu',iseed=12)
        ! 2)  Leave-n-out, where n=0.1*ntot
        CALL leave_some_out(0.1,'DATA_ws_Sc.dat','gp_ws_Sc_J26.emu','lnu_ws_Sc_J26.dat',iseed=12)
        ! 3)  IFS predictions
        CALL IFS_sample_full(output='test_emu_sample_all_Sc_J26.dat')
    ELSEIF (.FALSE.) THEN
        ! Combined night and day data from Nordling et al. (2024) - compile with "n24" flag
        ! 1) Training
        CALL test_emulator('DATA_merge_w',' ','test_wp_N24.dat','gp_wp_N24.emu',iseed=12)
        ! 2)  Leave-n-out, where n=0.1*ntot
        CALL leave_some_out(0.1,'DATA_merge_w','gp_wp_N24.emu','lnu_wp_N24.dat',iseed=12)
        ! 3)  IFS predictions
        CALL IFS_sample_full(output='test_emu_sample_all_N24.dat')
    ELSEIF (.FALSE.) THEN
        ! Training data by Noora Hyttinen from 2026 - compile with "h26" flag
        ! 1) Training
        CALL test_emulator('DATA_Noora_ws.dat',' ','test_ws_H26.dat','gp_ws_H26.emu',iseed=12)
        ! 2)  Leave-n-out, where n=0.1*ntot
        CALL leave_some_out(0.1,'DATA_Noora_ws.dat','gp_ws_H26.emu','lnu_ws_H26.dat',iseed=12)
        ! 3)  IFS predictions
        CALL IFS_sample_full(output='test_emu_sample_all_H26.dat')
    ENDIF
END PROGRAM me_main
