#!/usr/bin/env groovy

//inspired by: https://www.christopherrung.com/2017/05/04/slack-build-notifications/

import groovy.json.JsonOutput

def slackNotificationChannel = "dispatch-commits"
def commit = ""
def author = ""
def modifiedFiles = ""

def notifySlack(text, channel, attachments) {

    def slackURL = sh(returnStdout: true, script: 'less /mn/stornext/d9/mikolajs/Jenkins_CI/webhook_dispatch.cf').trim()

    def payload = JsonOutput.toJson([text: text,
        channel: channel,
        username: "Jenkins",
        attachments: attachments
    ])

    sh "curl -X POST --data-urlencode \'payload=${payload}\' ${slackURL}"
}

def uploadSummary(channel,file_name) {

   def botToken = sh(returnStdout: true, script: 'less /mn/stornext/d9/mikolajs/Jenkins_CI/app_bot_token.cf').trim()

   sh "curl -F file=@${file_name} -F \"initial_comment=Error Details\" -F channels=${channel} -H \"Authorization: Bearer ${botToken}\" https://slack.com/api/files.upload"

}

def getGitAuthor = {
    commit = sh(returnStdout: true, script: 'git rev-parse --short HEAD')
    author = sh(returnStdout: true, script: "git --no-pager show -s --format='%an' ${commit}").trim()
    author = "`"+author+"`"
    commit = "`"+commit.trim()+"`"
}

def getModifiedFiles = {
    modifiedFiles = sh(returnStdout: true, script: 'git diff --name-only HEAD HEAD~1 | head -n 10').trim() 
}

node {

       stage('Checkout'){
          checkout scm
       }

    try{
 
       stage("Post to Slack") {

           getGitAuthor()
           getModifiedFiles()
           
           def jobName = env.JOB_NAME.replace("/"," ") //prepare variable for use with sed command below 
               
           //load environment, create file with interesting information and update CI repo for website
           sh """#!/bin/bash -el
           python bifrost_test.py --setup_file experiments/tests_setup.yml &> output_file.txt

           """
 
           def testSummary = sh(returnStdout: true, script:'less summary.txt').trim()        

           testSummary = testSummary.replace("'","") //remove single quotes for now
           testSummary = "```" + testSummary + "```"

           notifySlack("", slackNotificationChannel, [
              [
                 title: "${env.JOB_NAME}, build #${env.BUILD_NUMBER}",
                 text: "Dear ${author}, congratulations your commit ${commit} passes or fails gracefully :+1:",
                 color: "#0339fc",
                 fields:[
                    [ 
                      title: "Modified Files",
                      value: "${modifiedFiles}",
                      short: false
                    ],
                    [
                       title: "Test Results",
                       value: "${testSummary}",
                       short: false
                    ]
                 ]
              ]
           ])

       }
    } catch(e){
 
        sh 'echo error!'
     
  
        def testSummary = sh(returnStdout: true, script:'less summary.txt').trim()    
     
        testSummary = testSummary.replace("'","") //remove single quotes for now
        testSummary = "```" + testSummary + "```"

        notifySlack("", slackNotificationChannel, [
           [
              title: "${env.JOB_NAME}, build #${env.BUILD_NUMBER}",
              text: "Hello ${author} you might want to take a closer look to ${commit} ",
              color: "#ffd700",
              fields:[
                 
                 [
                  title: "Jenkins console output",
                  value: "<${env.BUILD_URL}|Open>",
                  short: true
                 ],
                 [ 
                   title: "Modified Files",
                   value: "${modifiedFiles}",
                   short: false
                 ],

                 [
                    title: "Test Results",
                    value: "${testSummary}",
                    short: false
                 ]
              ]
           ]
        ])
       
        uploadSummary(slackNotificationChannel,'output_file.txt') 
    }
}


