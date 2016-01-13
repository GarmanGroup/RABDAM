import cProfile
cProfile.run('cambda("2BN3")', 'cambda.profile')
import pstats
stats = pstats.Stats('cambda.profile')
stats.strip_dirs().sort_stats('time').print_stats()
