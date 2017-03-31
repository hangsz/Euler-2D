
program main
    use Print_data
    implicit none

    call Grid        !网格读取，并求取需用几何参数
    call Flow_init   !流场初始化
    call Solver
    
    call Cp_cal
    call Output
end program
