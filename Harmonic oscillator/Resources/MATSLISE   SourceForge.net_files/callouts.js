/*global jQuery */

jQuery(function($) {
    function send_dismiss() { $.ajax({url: 'dismiss_callout' }); }

    $('a#edit-screenshots.add-button').qtip({
        content: {
            text:$('#screenshots-callout-content'),
            button:true
        },
        show: {
            ready: true,
            event: false
        },
        hide: {
            event: 'click',
            fixed: true
        },
        position: {
            my: 'bottom left',
            at: 'top right'
        },
        style: {
            classes: 'qtip-blue qtip-shadow'
        },
        events: {
            hide: send_dismiss
        }
    });

    $('a#edit-features.add-button').qtip({
        content: {
            text:$('#features-callout-content'),
            button:true
        },
        show: {
            ready: true,
            event: false
        },
        hide: {
            event: 'click',
            fixed: true
        },
        position: {
            my: 'bottom left',
            at: 'top right'
        },
        style: {
            classes: 'qtip-blue qtip-shadow'
        },
        events: {
            hide: send_dismiss
        }

    });

    $('.sfdl, .sfdl small, .info-circle').qtip({
    });
});

