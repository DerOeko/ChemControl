# Motivational Cognitive Control – ChemControl Models

In order to investigate how the controllability of stressful outcomes affect the recruitment of dopaminergic and serotonergic pathways using computational modeling, we built several Rescorla-Wagner models to fit behavioral data retrieved from a study manipulating controllability.

In total, we've build 4 RW models based on Guitart—Masip et al. 2012, and implemented several models of our own to model a separate controllability estimate.

## Using this GitHub repository:

* Clone this GitHub repository.
* Create a "github_config.m" file that creates a folderPath string that points to your data. E.g. like this:

"% Local path on your machine
folderPath = "/project/path/to/your/data/";"
* Specify hyperparameters and metaparameters in the Config class. This class is passed to all models and centralizes all necessary parameters.
* Done!

## Contact
For questions, contact me, at samuelgerrit.nellessen@gmail.com.
