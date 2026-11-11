# Accept only checkpoints from the verified pre-separation implementation.
null_checkpoint_legacy <- function(old, current) {
  previous <- c('8c9834479d59d2ecf801c2acf1e6fee0', 'd9d70da413655505a6820497c64020c6',
                '46ea4ebc2957f0260452f194ad2edc95', 'c025a0e4909f1cb439d273d9724151cc',
                'e79a6ee06e7836a2dacf2afd1b22e6c2')
  old_separated <- c('8c9834479d59d2ecf801c2acf1e6fee0', 'd9d70da413655505a6820497c64020c6',
                 '46ea4ebc2957f0260452f194ad2edc95', 'dfc78de1b3b6eb84e5e2e1ae4543fddd',
                 '0dbe49c4633b1dfb8b2f75043fdc152c')
  separated <- old_separated
  separated[5] <- 'cdcaac7f7da81730c1f4b92568006147'
  same_scripts <- (identical(old$scripts, previous) || identical(old$scripts, old_separated)) &&
    identical(current$scripts, separated)
  old$scripts <- current$scripts
  same_scripts && identical(old, current)
}
